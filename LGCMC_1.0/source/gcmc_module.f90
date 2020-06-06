!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             This file is a part of the code            !
!                ╦    ╔═╗  ╔═╗  ╔╦╗  ╔═╗                 !
!                ║    ║ ╦  ║    ║║║  ║                   !
!                ╩═╝  ╚═╝  ╚═╝  ╩ ╩  ╚═╝                 !
!       A Lattice Grand Canonical Monte Carlo code       !
!                Written by Arpan Kundu                  !
!                        at                              !
!            Humboldt University of Berlin               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 MODULE gcmc_module
 USE variable
 USE site_finder_module, ONLY : site_finder
 USE energy_calculation, ONLY : energy
 USE random_number_generator, ONLY : initialize_rng,ran_num
      
!      USE calc_ads_en, ONLY : ads_en

 DOUBLE PRECISION :: kt   !thermal energy
 DOUBLE PRECISION :: p    !pressure
 DOUBLE PRECISION,ALLOCATABLE::phi(:),f(:),x(:) !fugacity coeff, fugacity, composition
 DOUBLE PRECISION, PRIVATE :: eps = 1.0E-16     !smallest number
 INTEGER :: max_ex_sel_trial=10                 !maximum number of allowed trial to select an exchange process
 LOGICAL :: virtual_count=.FALSE.               !to initiate when the virual_states counting is needed.
 INTEGER :: nvirtual=0                            !number of virtual sites at a particular moment
 INTEGER :: i_kw                                ! number of kawasaki dynamics counter


 CONTAINS

 SUBROUTINE gcmc(avgocc,dist_f,corr_f,avgocc_evn_odd)
 USE variable

 IMPLICIT NONE
! ### Declearation of variables and arrays:
 DOUBLE PRECISION avgocc(ncomp,nsite_type)
 DOUBLE PRECISION avgocc_evn_odd(ncomp,nsite_type,4) ! (X=Even,Y=Even)=>1,(Even,Odd)=>2,(Odd,Even)=>3,(Odd,Odd)=>4
 DOUBLE PRECISION dist_f(ncomp,nsite_type,corx1:corx2,cory1:cory2,ncomp,nsite_type,9) !X=Even =>5, Y=Even =>6, X=Odd =>7, Y=Odd =>8 Total => 9

 DOUBLE PRECISION corr_f(ncomp,nsite_type,corx1:corx2,cory1:cory2,ncomp,nsite_type,9)

 INTEGER imcstep,exchange_attempt,icomp,ixsite,iysite,isite_type !counters
 INTEGER isite_type1,isite_type2,ix,iy,icomp1,icomp2             !counters      
 INTEGER evn_odd,evn_odd2                                        ! (X=Even,Y=Even)=>1,(Even,Odd)=>2,(Odd,Even)=>3,(Odd,Odd)=>4
 INTEGER jxsite,jysite                                           !counter
 INTEGER disp_attempt
 INTEGER config_taken                                             !number of configurations taken in sum
 INTEGER tempocc(0:ncomp,1:nsite_type)
 INTEGER(KIND=8) sumocc(ncomp,nsite_type)
 INTEGER(KIND=8) sum_corr(1:ncomp,1:nsite_type,corx1:corx2,cory1:cory2,1:ncomp,1:nsite_type,4)
 DOUBLE PRECISION old_energy, new_energy

 INTEGER(KIND=8) :: n_tot_kw,n_tot_kw_non_dg,n_acc_kw            !Total no all, non degenerate Kw and no of accepted Kw move
 INTEGER(KIND=8) :: n_tot_gl,n_acc_gl                            ! Total no and no of accepted Glauber move
 DOUBLE PRECISION :: p_acc_kw, p_acc_gl,avg_kw_non_dg,avg_kw_real !Glauber and Kawasaki acceptance, average non-degenerate kw per MC cycle
 DOUBLE PRECISION :: avg_kw                                       !Average effective diffusion tried
 LOGICAL :: disp_possible                 !If all sites have same configuration then disp is not possible

 INTEGER :: s_v                           ! survived virtual state
 INTEGER(KIND=8) :: n_kw_real                     !  real kw attempted
 INTEGER :: kw_degenerate                 ! no of degenerate trial tried per mcstep

 LOGICAL :: acc_gl,acc_kw,site_selection      ! accept/reject of Glauber and Kawasaki moves.
 INTEGER :: eqbm_tracker, prod_tracker        ! Progress tracker


 eqbm_tracker=neqbmstep/20
 prod_tracker=nsimulstep/20

 s_v=0

!****  open(203,file='last_config',status='unknown')
!****  write(203,*) 'Configuration after',ntotalstep
!****  write(203,'(1X,A5,2x,A5,2x,A9,2x,6A)') 'XSITE','YSITE','SITE_TYPE',



 kt=8.31*t/1000.0d0

 CALL setup

 config_taken=0

 n_tot_kw = 0
 n_tot_kw_non_dg=0
 n_acc_kw = 0
 n_acc_gl = 0
 n_kw_real = 0
 n_tot_gl=ntotalstep*nsite*nsite_type

!**************** Test interaction print*************************
!*       do icomp1=1,ncomp
!*         do isite_type1=1,nsite_type
!*           do ix=-1*xcut,xcut
!*              do iy=-1*ycut,ycut
!*                do icomp2=1,ncomp
!*                  do isite_type2=1,nsite_type
!*
!*       write(*,*) icomp1,icomp2,isite_type1,isite_type2,ix,iy,
!*     x j(icomp1,isite_type1,ix,iy,icomp2,isite_type2),
!*     x jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)
!*                   enddo
!*               enddo
!*             enddo
!*            enddo
!*            enddo
!*            enddo
!******************************************************************

 do icomp1=1,ncomp
    do isite_type1=1,nsite_type
       sumocc(icomp1,isite_type1)=0
       do ix=corx1,corx2
          do iy=cory1,cory2
             do icomp2=1,ncomp
                do isite_type2=1,nsite_type
                   do evn_odd=1,4

                      sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)=0

                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo


 do 100 imcstep=1,ntotalstep
!****   write(*,*) 'MC step running',imcstep
  exchange_attempt=0
  tempocc(:,:) = 0
  do while (exchange_attempt.lt.(nsite*nsite_type))
     call exchange(acc_gl)
     exchange_attempt=exchange_attempt+1
     if(acc_gl) then
       n_acc_gl=n_acc_gl+1
     endif
  enddo

  if((imcstep.eq.disp_start).and.(constraint_file_needed)) then
    call make_virtual_list
  endif

!*****  write(*,*) 'MC step running', imcstep
 
  i_kw = 0                                            !Setting Kawasaki dynamics counter to 0
  disp_possible=.TRUE.


1200 continue

!****   write(*,*) 'MCstep/Disp',imcstep,(displacement.and.disp_possible)

  if(displacement.and.disp_possible) then
    if(imcstep.ge.disp_start) then
!*****        write(*,*) 'Starting displacement @',imcstep
!****        Do disp_attempt=1,maxdisp
       kw_degenerate=0
       do while (.not.((nvirtual.eq.0).and.(i_kw.ge.maxdisp)))
!            Do while (i_kw.lt.maxdisp)
!               write(*,'(I6)',advance='no') i_kw
          call translate(acc_kw,site_selection)
          if(acc_kw) then
            n_acc_kw=n_acc_kw+1
          endif
          if(site_selection) then
            n_kw_real=n_kw_real+1
          else
            kw_degenerate=kw_degenerate+1
          endif  
       Enddo
!*     write(*,*) "Attempted Kw=",i_kw,imcstep,nvirtual
        n_tot_kw=n_tot_kw+int8(i_kw)
        n_tot_kw_non_dg = n_tot_kw_non_dg+int8(i_kw)-int8(kw_degenerate)
    endif
  endif

  !Carriage control: Reporting progress bar for equilibration/simulation.
  open (unit=6, carriagecontrol='fortran')
  if (imcstep.eq.1) then
       write(*,'(A)',advance='no') " Equilibration Running: |"

  else if (imcstep.le.neqbmstep) then
       if (Mod(imcstep,eqbm_tracker).eq.0) then
            write(*,'(A)',advance='no') "+"      
       endif
       if (imcstep.eq.neqbmstep) then
            write(*,'(A)') "|"
            write(*,'(A)',advance='no') " Simulation Running:    |"
       end if

  else if ((imcstep.gt.neqbmstep).and.(imcstep.le.ntotalstep)) then
       if (Mod(imcstep-neqbmstep,prod_tracker).eq.0) then
            write(*,'(A)',advance='no') "+"
       endif
       if (imcstep.eq.ntotalstep) then
          write(*,'(A)') "|"
       endif           
 endif
! Writing the statistics at end
  if(imcstep.eq.ntotalstep) then
    p_acc_gl=dble(n_acc_gl)/dble(n_tot_gl)
    write(*,*) "Exchange with reservoir acceptance ratio",p_acc_gl

    if(displacement) then
    p_acc_kw=dble(n_acc_kw)/dble(n_tot_kw_non_dg)
    avg_kw_non_dg=dble(n_tot_kw_non_dg)/dble(ntotalstep-disp_start+1)
    avg_kw=dble(n_tot_kw)/dble(ntotalstep-disp_start+1)
    avg_kw_real = dble(n_kw_real)/dble(ntotalstep-disp_start+1)

   write(*,*) "Exchange with reservoir acceptance ratio",p_acc_gl
   write(*,*) "Average number of effective diffusion attempted",NINT(avg_kw)
   write(*,*) "Average number of real diffusion attempted",NINT(avg_kw_real)
   write(*,*) "Average number of non-degenerate effective diffusion attempted", NINT(avg_kw_non_dg)
   write(*,*) "Non-degenerate diffusion acceptance ratio",p_acc_kw
   endif
  endif

  if (imcstep.gt.neqbmstep) then
    if(Mod(imcstep-neqbmstep,tdump).eq.0) then

      if(nvirtual.gt.0) then
!*             write(*,*) imcstep,i_kw ! nvirtual
        s_v=s_v+nvirtual
      endif

      config_taken=config_taken+1
      do ixsite=1,nxsite
        do iysite=1,nysite
          do icomp1=1,ncomp
            do isite_type1=1,nsite_type
              sumocc(icomp1,isite_type1)=sumocc(icomp1,isite_type1)+ kronecker_delta(icomp1,confg(ixsite,iysite,isite_type1))
          
               if(.not.corr_func) then
                 cycle
               endif

! Calculation of Correlation function and other statistics

                if((Mod(ixsite,2).eq.0).and.(Mod(iysite,2).eq.0)) then
                      evn_odd=1
                else if ((Mod(ixsite,2).eq.0).and.(Mod(iysite,2).eq.1)) then
                      evn_odd=2
                else if ((Mod(ixsite,2).eq.1).and.(Mod(iysite,2).eq.0)) then
                      evn_odd=3
                else
                      evn_odd=4
                endif

                  do ix=corx1,corx2
                     do iy=cory1,cory2

                        call site_finder(ixsite,iysite,ix,iy,jxsite,jysite)
                        do icomp2=1,ncomp
                           do isite_type2=1,nsite_type
                          
                              sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= & 
                              & sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd) + & 
                              & kronecker_delta(icomp1,confg(ixsite,iysite,isite_type1))* &
                              & kronecker_delta(icomp2,confg(jxsite,jysite,isite_type2))
                                         
                           enddo
                        enddo
                     enddo
                  enddo
            enddo 
          enddo
        enddo
      enddo
    endif
  endif

100   continue

 do icomp1=1,ncomp
    do isite_type1=1,nsite_type
       avgocc(icomp1,isite_type1)=dble(sumocc(icomp1,isite_type1))/dble(nsite*config_taken)
    enddo
 enddo

 
 if(corr_func) then
 do icomp1=1,ncomp
  do isite_type1=1,nsite_type
   do evn_odd =1,9

      if(evn_odd.le.4) then
      avgocc_evn_odd(icomp1,isite_type1,evn_odd)= &
 &    sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,evn_odd)/dble(nsite*config_taken)
      endif     
 
    do ix=corx1,corx2
     do iy=cory1,cory2
      do icomp2=1,ncomp
       do isite_type2=1,nsite_type

      if (evn_odd.le.4) then
      dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= & ! P(B|A) where B=icomp2,isite_type2 and A=icomp1,isite_type1 
   &  dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd))/ & 
   &  dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,evn_odd))

      
      else if (evn_odd.eq.5) then                                  ! X=Even
      dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &
  &   dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,1)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,2))/ &
  &   dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,1)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,2))

      else if (evn_odd.eq.6) then                                  ! Y=Even
      dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &
  &   dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,1)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,3))/ &
  &   dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,1)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,3))

      else if (evn_odd.eq.7) then                                  ! X = Odd
      dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &
  &   dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,3)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,4))/ &
  &   dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,3)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,4))

      else if (evn_odd.eq.8) then                                  ! Y = Odd
      dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &
  &   dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,2)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,4))/ &
  &   dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,2)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,4))

    else
       dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &
  &    dble(sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,1)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,2) &
    &       +sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,3)+sum_corr(icomp1,isite_type1,ix,iy,icomp2,isite_type2,4))/ &
  &    dble(sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,1)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,2) &
    &       +sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,3)+sum_corr(icomp1,isite_type1,0,0,icomp1,isite_type1,4))

    endif

   
      corr_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)= &  !P(B|A)/P(B)=P(A.B)/P(A).P(B)  is the correlation function
   &    (dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,evn_odd)/avgocc(icomp2,isite_type2))-1
      enddo
       
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif

 if(constraint_file_needed) then
 write(*,*) 'survived virtual state', s_v
 endif
!****       do icomp=1,ncomp
!****          do isite_type=1,nsite_type
!****             write(*,*) avgocc(icomp,isite_type)
!****          enddo
!****       enddo
 RETURN
 END SUBROUTINE GCMC



 SUBROUTINE setup
 IMPLICIT NONE
 INTEGER ixsite,iysite,isite_type,icomp !counter variables

 If(.not.allocated(confg)) then
   ALLOCATE (confg(1:nxsite,1:nysite,1:nsite_type))
 Endif

 call initialize_rng(iseed)

 do icomp=1,ncomp
    f(icomp)=phi(icomp)*x(icomp)*p
 end do
    f(0)=1

 do ixsite=1,nxsite
    do iysite=1,nysite
       do isite_type=1,nsite_type
           confg(ixsite,iysite,isite_type)=0  !initialization of all confg at 0....clean surface
       enddo
    enddo
 enddo

 RETURN
 END SUBROUTINE setup

 ! Subroutine to exchange particle with reservoir
 SUBROUTINE exchange(accept)
 IMPLICIT NONE
 INTEGER ixwork,iywork,iswork  !working site index (x,y,site_type)
 INTEGER old_confg,new_confg
 DOUBLE PRECISION rand1,rand2,rand3,rand4,rand_site_type
 DOUBLE PRECISION dele, prob_acc
 LOGICAL :: accept

! Variable to bookkeep constraint sites
 TYPE(DEPENDENT) :: vir_loc1(2*max_dependent)
 TYPE(DEPENDENT) :: vir_cre(2*max_dependent)
 TYPE(DEPENDENT) :: vir_anh(2*max_dependent)
 INTEGER :: nv1, n_cre, n_anh

 accept=.FALSE.

2001  call ran_num(rand1)
 ixwork=int(rand1*dble(nxsite))+1
 if((ixwork.lt.1).or.(ixwork.gt.nxsite)) then
!****        write(*,*) "ixwork out of bound in exchange",ixwork
!****        write(*,*) "rand1=",rand1
   go to 2001
 endif

2002  call ran_num(rand2)
 iywork=int(rand2*dble(nysite))+1
 if((iywork.lt.1).or.(iywork.gt.nysite)) then
!****        write(*,*) "iywork out of bound in exchange",iywork
!****        write(*,*) "rand2=",rand2
   go to 2002
 endif

 if(nsite_type.gt.1) then
2003    call ran_num(rand_site_type)
   iswork=int(rand_site_type*dble(nsite_type))+1
   if((iswork.lt.1).or.(iswork.gt.nsite_type)) then
!****           write(*,*) "iswork out of bound in exchange",iswork
!****           write(*,*) "rand=",rand_site_type
      go to 2003
   endif
 else
   iswork=1
 endif

 old_confg=confg(ixwork,iywork,iswork)

! Selecting a new configuration (except the old one) at random.      
2004  call ran_num(rand3)
 new_confg=int(rand3*dble(ncomp))+1
 if((new_confg.lt.1).or.(new_confg.gt.ncomp)) then
!****          write(*,*) "new_confg out of bound in exchange",new_confg
!****          write(*,*) "rand3=",rand3
 go to 2004
 endif

 if(new_confg.eq.old_confg) then
   new_confg=0
 endif
 dele=energy(ixwork,iywork,iswork,new_confg)-energy(ixwork,iywork,iswork,old_confg)
 
 call ran_num(rand4)

 prob_acc=exp(-1.0*dele/kt)*f(new_confg)/f(old_confg)

 if (rand4.le.prob_acc) then
    if(virtual_count) then
     call search_virtual_local(ixwork,iywork,iswork,nv1,vir_loc1) !Making local list before move
    endif

    confg(ixwork,iywork,iswork)=new_confg                          !Exchange done
    
    if(virtual_count) then 
     call virtual_status_change(ixwork,iywork,iswork,vir_loc1,nv1,vir_cre, n_cre,vir_anh,n_anh) !Calculating creation and annihilation 
     call update_virtual_list(n_cre,vir_cre,n_anh,vir_anh)       !updating the list of virtual sites

!**************************DEBUG******************************************************************
!****         if((n_anh.gt.0).or.(n_cre.gt.0)) then
!****              write(*,*)'Gl', n_cre,n_anh,nvirtual
!****              write(*,*)'Gl',ixwork,iywork,iswork,old_confg,'->',
!****     x        confg(ixwork,iywork,iswork)
!****         endif
!*************************************************************************************************

    endif

    accept = .TRUE.
 endif


 RETURN
 END SUBROUTINE exchange


 SUBROUTINE translate(accept,site_selection)
 IMPLICIT NONE

 TYPE KEX                        !Kawasaki Exchange
 SEQUENCE
     INTEGER :: x,y,st            !xposition,yposition,site_type
     INTEGER :: init,final        !initial and final configuration
     DOUBLE PRECISION :: bf       ! Boltzman factor
     LOGICAL :: enable            ! if 'T' process is enabled else disabled
 ENDTYPE KEX

 
 INTEGER :: i,j                   !general counter
 INTEGER :: ixwork1,iywork1,ixwork2,iywork2,iswork1,iswork2
 INTEGER :: ivirtual
 INTEGER :: ix2,iy2,is2,ex_sel_trial
 INTEGER :: ixw,iyw                    ! site location needed for calculation of W(o) -> rosenbluth factor
 INTEGER :: ex_index,stor_ind              !exchange index inside the displacement box and storage index
 INTEGER :: config1,config2
 INTEGER :: n_dis,n_en                     !number of disabled process and enabled process
 INTEGER :: nsite_choice1,nsite_choice2      !number of possibility of chosing first site
 INTEGER :: nvirtual_temp                    !temporary nvirtual
 INTEGER :: chosen_process             !chosen_exchange_process
 DOUBLE PRECISION :: rand1,rand2,rand3,rand4,rand5,rand6
 DOUBLE PRECISION :: rand_site1
 DOUBLE PRECISION :: e_old,e_new,dele,prob_acc,prob_acc2
 DOUBLE PRECISION :: cbf(0:disp_box)       !cumulative boltzman factor
 DOUBLE PRECISION :: sum_bf              ! sum of Boltzman factors
 DOUBLE PRECISION :: wn,wo               !Rosenbluth factors
 LOGICAL :: site_selection,accept
 DOUBLE PRECISION p_corr,p_reject        !correction term of probability and probability of rejecting a particular site
 INTEGER :: inc_kw                       !increment of Kawasaki exchange

 TYPE (KEX) :: ex_info(disp_box) 
 TYPE (KEX) :: ex_tmp                    !temporary to store data

 ! Variable to bookkeep constraint sites
 TYPE(DEPENDENT), DIMENSION(2*max_dependent) :: vir_loc1, vir_loc2 
 TYPE(DEPENDENT), DIMENSION(2*max_dependent) :: vir_cre1,vir_cre2
 TYPE(DEPENDENT), DIMENSION(2*max_dependent) :: vir_anh1,vir_anh2
 TYPE(DEPENDENT), DIMENSION(2*max_dependent) :: vir_cre,vir_anh    !Total virtual sites created or annihilated in totality
 INTEGER :: nv1,nv2,n_cre1,n_cre2,n_cre,n_anh1,n_anh2,n_anh


 accept=.FALSE.
 site_selection=.FALSE.
 ex_sel_trial=0
 
 

!      Do While (.NOT.site_selection)                    !We do not use the do loop anymore as failure of site selection will be treated as degenerate diffusion/ degenerate trial.

!       if(ex_sel_trial.ge.max_ex_sel_trial) then
!         exit
!       endif

!       ex_sel_trial=ex_sel_trial+1

if(nvirtual.gt.0) then                               !if there are virtual sites present then always a virtual site be selected
  nsite_choice1=nvirtual                              !as first site and the corresponding acceptance probabilities will be modified
                                                       !due to the bias introduced
2008  call ran_num(rand1)         
  ivirtual=int(rand1*dble(nvirtual))+1
  if(ivirtual.gt.nvirtual) then
    goto 2008
  endif

  ixwork1 = virtual_sites(ivirtual)%x
  iywork1 = virtual_sites(ivirtual)%y
  iswork1 = virtual_sites(ivirtual)%site

!****************DEBUG**********************************
!******         write(*,*) ixwork1,iywork1,iswork1
!*******************************************************   
else
  nsite_choice1=nsite*nsite_type  
 
2005  call ran_num(rand1)
  ixwork1=int(rand1*dble(nxsite))+1
  if((ixwork1.lt.1).or.(ixwork1.gt.nxsite)) then
!****           write(*,*) "ixwork1 out of bound in translate",ixwork1
!****           write(*,*) "rand1=",rand1
  go to 2005
  endif

2006 call ran_num(rand2)
  iywork1=int(rand2*dble(nysite))+1
  if((iywork1.lt.1).or.(iywork1.gt.nysite)) then
!****           write(*,*) "iywork1 out of bound in translate",iywork1
!****           write(*,*) "rand2=",rand2
    go to 2006
  endif

  if(nsite_type.gt.1) then
2007 call ran_num(rand_site1)
     iswork1=int(rand_site1*dble(nsite_type))+1
     if((iswork1.lt.1).or.(iswork1.gt.nsite_type)) then
!****              write(*,*) "iswork1 out of bound in translate",iswork1
!****              write(*,*) "rand=",rand_site1
        go to 2007
     endif
  else
      iswork1=1
  endif
 
endif 
          
!Scanning the all displacements in the box

   ex_index=0          !scanning each exchange process within the box
   n_dis = 0           ! Number of disabled exchange
   wn = 0.0d0          ! initialize W(n): rosenbluth factor for flow towards new state
   
   do ix2=-1*x_disp,x_disp
     do iy2=-1*y_disp,y_disp
       do iswork2=1,nsite_type
          if((ix2.eq.0).and.(iy2.eq.0).and.(iswork2.eq.iswork1)) then   !Removing self exchange :P
             cycle
          endif

          call site_finder(ixwork1,iywork1,ix2,iy2,ixwork2,iywork2)

          ex_index=ex_index+1

          config1=confg(ixwork1,iywork1,iswork1)
          config2=confg(ixwork2,iywork2,iswork2)
          
          if (config1.eq.config2) then 
             n_dis=n_dis+1
             stor_ind=disp_box-n_dis+1                            !Disabled processes are stored in the array from back side
             ex_info(stor_ind)%enable = .FALSE.
             ex_info(stor_ind)%bf = 0.0d0   !1.0d0
             go to 2100
          endif

          e_old=energy(ixwork1,iywork1,iswork1,config1)+energy(ixwork2,iywork2,iswork2,config2)
          e_new=energy(ixwork1,iywork1,iswork1,config2)+energy(ixwork2,iywork2,iswork2,config1)
          dele=e_new-e_old 


          if(((e_old.gt.h_cons_min).and.(e_new.gt.h_cons_min)).or.(abs(dele).lt.eps)) then
             n_dis=n_dis+1
             stor_ind=disp_box-n_dis+1                           !Disabled processes are stored from the back side of the array
             ex_info(stor_ind)%enable = .FALSE.
             ex_info(ex_index)%bf = 0.0d0   !1.0d0
             goto 2100
          endif 

          stor_ind=ex_index-n_dis                                !Enabled processes are stored in forward direction
          ex_info(stor_ind)%bf=exp(-1.0*dele/kt)
          ex_info(stor_ind)%enable = .TRUE.

2100      continue

          ex_info(stor_ind)%x  = ixwork2
          ex_info(stor_ind)%y  = iywork2
          ex_info(stor_ind)%st = iswork2
          ex_info(stor_ind)%init  = config2
          ex_info(stor_ind)%final = config1
          
          wn=wn+ex_info(stor_ind)%bf                  !calculation of rosenbluth factor at end

       enddo
     enddo
   enddo  

   if (n_dis.eq.disp_box) then
!*            write(*,*) 'no possible displacement', nvirtual
       i_kw = i_kw+1                    ! One attempted trial
       RETURN
!*           cycle
   endif

!sorting the possible processes

   n_en=disp_box-n_dis                                !number of enabled process

   do i=1,n_en-1
     do j=i+1,n_en
        if(ex_info(i)%bf.lt.ex_info(j)%bf) then
          ex_tmp=ex_info(i)
          ex_info(i)=ex_info(j)
          ex_info(j)=ex_tmp
        endif
     enddo
   enddo

!calculating cumulative probability
  sum_bf=0.0d0
  cbf(0)=0.0d0
  do i=1,n_en
    sum_bf=sum_bf+ex_info(i)%bf
    cbf(i)=sum_bf
  enddo            

!selecting a site
  call ran_num(rand3)
  rand4=rand3*sum_bf  

  i=1
  if(n_en.gt.1) then
   do while(.not.((rand4.le.cbf(i)).and.(rand4.gt.cbf(i-1)))) 
     i=i+1
     if(i.eq.n_en) then
!*            write(*,*) i
       exit
     endif
   enddo
  endif   

  chosen_process = i     
  site_selection=.TRUE.
  ixwork2=ex_info(chosen_process)%x
  iywork2=ex_info(chosen_process)%y
  iswork2=ex_info(chosen_process)%st


  config1=ex_info(chosen_process)%final
!Calculating Rosenbluth factor for the cell arround chosen site (wo)

!************ Debug**********************
!*       write(*,*) 'process_list'
!*       write(*,*) (ex_info(i)%bf,i=1,n_en)
!*       write(*,*) 'chosen process'
!*       write(*,*) ixwork1,iywork1,iswork1,ixwork2,iywork2,iswork2
!*       write(*,*) ex_info(chosen_process)%final, 
!*     x            ex_info(chosen_process)%init
!*       write(*,*) chosen_process
!*       write(*,*) cbf(chosen_process)/sum_bf
!*       write(*,*) ex_info(chosen_process)%bf
!       
!************ Debug************************

  wo=0.0d0
  do ix2=-1*x_disp,x_disp
   do iy2=-1*y_disp,y_disp
     do is2=1,nsite_type
          if((ix2.eq.0).and.(iy2.eq.0).and.(is2.eq.iswork2)) then
             cycle
          endif 
          call site_finder(ixwork2,iywork2,ix2,iy2,ixw,iyw)
          if((ixw.eq.ixwork1).and.(iyw.eq.iywork1)) then
            wo=wo+(1/ex_info(chosen_process)%bf)
            cycle
          endif

          config2=confg(ixw,iyw,is2)

          if(config1.eq.config2) then
            wo=wo+0.0d0   !1.0d0
            cycle
          endif

          e_old=energy(ixwork2,iywork2,iswork2,config1)+energy(ixw,iyw,is2,config2)
          e_new=energy(ixwork2,iywork2,iswork2,config2)+energy(ixw,iyw,is2,config1)
          dele=e_new-e_old

          if(((e_old.gt.h_cons_min).and.(e_new.gt.h_cons_min)).or.(abs(dele).lt.eps)) then
            wo=wo+0.0d0      !1.0d0
            cycle
          endif

          wo=wo+exp(-1.0d0*dele/kt)

       enddo
     enddo
  enddo  
          
!*       prob_acc=(wn*dble(nsite_choice))/
!*     x  (wo*ex_info(chosen_process)%bf*dble(disp_box))


  call ran_num(rand5)

  if(virtual_count) then

    call search_virtual_local(ixwork1,iywork1,iswork1,nv1,vir_loc1)!Making local list before move for site1
    call search_virtual_local(ixwork2,iywork2,iswork2,nv2,vir_loc2)!Making local list before move for site2

    confg(ixwork1,iywork1,iswork1)=ex_info(chosen_process)%init   !Kawasaki exchange is done here
    confg(ixwork2,iywork2,iswork2)=ex_info(chosen_process)%final  ! -do-

    call virtual_status_change(ixwork1,iywork1,iswork1,vir_loc1,nv1,vir_cre1, n_cre1,vir_anh1,n_anh1)  !Calculating creation and annihilation for site 1 
    call virtual_status_change(ixwork2,iywork2,iswork2,vir_loc2,nv2,vir_cre2, n_cre2,vir_anh2,n_anh2)  !Calculating creation and annihilation for site 1 nv2,vir_cre2, n_cre2,vir_anh2,n_anh2) 
    
    call virtual_kawasaki(n_cre1,vir_cre1,n_cre2,vir_cre2,n_anh1,vir_anh1,n_anh2,vir_anh2,n_cre,vir_cre,n_anh,vir_anh)    ! Calculating total creation and annihilation.
   
    nvirtual_temp=nvirtual-n_anh+n_cre

    if(nvirtual_temp.ne.0) then
       i_kw = i_kw+1
       nsite_choice2=nvirtual_temp
    else
       if(nvirtual.gt.0) then
        p_reject=dble(nsite*nsite_type-1)/dble(nsite*nsite_type)   !probability of not selecting a site
        call ran_num(rand6)
        inc_kw = int(log(rand6)/log(p_reject))+1                   !inc_kw is random number from geometric distribution.
        i_kw=i_kw+inc_kw

!*****************************************************************************************************************
! If you want to stop the calculation @ maxdisp then use this block and 3rd do while loop in displacement block
!******************************************************************************************************************
!             if(i_kw.gt.maxdisp) then
!               write(*,*) "Move finished before relaxation",nvirtual
!               i_kw=maxdisp
!          confg(ixwork1,iywork1,iswork1)=ex_info(chosen_process)%final !configurations are restored as the move is rejected with modified probability
!          confg(ixwork2,iywork2,iswork2)=ex_info(chosen_process)%init  ! -do-
!              accept=.FALSE.
!              RETURN
!             endif
!***********************************************************************************************************************            
        nsite_choice2=1                                            !nsite*nsite_type artificially set it to 1 to avoid low acceptance probability.
       else
       i_kw=i_kw+1
       nsite_choice2=nsite*nsite_type
       endif
!*            if(nvirtual.gt.0) then
!*            write(*,*) 'wn,wo,cbf',wn,wo,ex_info(chosen_process)%bf
!*            endif
    endif

    p_corr=dble(nsite_choice1)/dble(nsite_choice2)
    prob_acc=p_corr*wn/(wo*ex_info(chosen_process)%bf)
   
    if(rand5.le.prob_acc) then
!*           if(nvirtual.gt.0) then
!*             write(*,*) 'Relaxation accepted'
!*           endif

     call update_virtual_list(n_cre,vir_cre,n_anh,vir_anh)         ! Updating the virtual list as already accepted
     accept=.TRUE.
!*          write(*,*) 'Move accepted',n_cre, n_anh
    else
     confg(ixwork1,iywork1,iswork1)=ex_info(chosen_process)%final !configurations are restored as the move is rejected with modified probability
     confg(ixwork2,iywork2,iswork2)=ex_info(chosen_process)%init  ! -do- 
     accept=.FALSE.

!*            if(nvirtual.gt.0) then
!*              write(*,*) 'Relaxn rejected'
!*            endif
    endif
    
  else
    i_kw = i_kw+1
    prob_acc=wn/(wo*ex_info(chosen_process)%bf)
    if(rand5.le.prob_acc) then
      confg(ixwork1,iywork1,iswork1)=ex_info(chosen_process)%init  !Kawasaki exchange is done here
      confg(ixwork2,iywork2,iswork2)=ex_info(chosen_process)%final ! -do-
      accept=.TRUE.
    endif
    
  endif
       
!*       i_kw=i_kw+inc_kw
!*******************************DEBUG***************************
!*       if((n_anh.gt.0).or.(n_cre.gt.0)) then
!***       write(*,*)'Kw',n_cre,n_anh,nvirtual
!*       endif
!*        write(*,*) 'Kw',ixwork1,iywork1,iswork1,
!*     x confg(ixwork2,iywork2,iswork2),'->',ixwork2,iywork2,iswork2,
!*     x confg(ixwork1,iywork1,iswork1)
!**************************************************************
 
!**************************************************************

!      enddo

!*******************DEBUG************************************ 
!*        write(*,*) 'wn',wn
!*        write(*,*) 'wo', wo
!*        write(*,*) 'acc_prob=', prob_acc
!*        write(*,*) accept
!
!      print *, "site1",ixwork1,iywork1,iswork1 , "site2", 
!     x         ixwork2, iywork2,iswork2
!         write(*,*) nvirtual
!******************DEBUG*************************************
 RETURN
 END SUBROUTINE


 SUBROUTINE make_virtual_list
 USE  calc_ads_en  !, only : constraint_satisfied
 IMPLICIT NONE
 INTEGER :: ii,jj,kk,ll,mm     !counter variables
 INTEGER :: ivirtual         !counter for virtual sites

 ivirtual=0
 if(.not.allocated(virtual_sites)) then
   allocate(virtual_sites(nsite*nsite_type))
 endif

 do ii = 1,nxsite
  do jj = 1,nysite
   do kk = 1,nsite_type
     if(.not.constraint_satisfied(ii,jj,kk,confg(ii,jj,kk))) then
!*****            write(*,*) "constraint false"
       ivirtual=ivirtual+1 
       virtual_sites(ivirtual)%x=ii
       virtual_sites(ivirtual)%y=jj
       virtual_sites(ivirtual)%site=kk
     endif


   enddo
  enddo
 enddo

 virtual_count=.TRUE.
 nvirtual=ivirtual

!***************Debug************************
!*      write(*,*) 'Vrtual list created ', nvirtual
!***************Debug*************************
 RETURN
 END SUBROUTINE make_virtual_list

 SUBROUTINE search_virtual_local(ix,iy,is,nv,vir_loc)
 USE calc_ads_en, only : constraint_satisfied
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ix,iy,is
 INTEGER :: jx,jy,js                                      ! Other site counter
 INTEGER :: ii,jj,kk                                      ! Directions 
 INTEGER :: iv                                          ! counter of virtual sites
 INTEGER :: id                                            ! counter for scanning dependent sites
 INTEGER :: nv                                            !number of virtual sites
 TYPE(DEPENDENT) :: vir_loc(2*max_dependent)     ! virtual_sites in local environment
   
 iv=0
 if(.not.constraint_satisfied(ix,iy,is,confg(ix,iy,is))) then
   iv=iv+1
   vir_loc(iv)%x=ix
   vir_loc(iv)%y=iy
   vir_loc(iv)%site=is
 endif

 if(no_of_dependent(is).gt.0) then
  do id = 1,no_of_dependent(is)

     ii = dependent_list(is,id)%x
     jj = dependent_list(is,id)%y
     js = dependent_list(is,id)%site

     call site_finder(ix,iy,ii,jj,jx,jy)
     if(.not.constraint_satisfied(jx,jy,js,confg(jx,jy,js))) then
       iv=iv+1           
       vir_loc(iv)%x=jx
       vir_loc(iv)%y=jy
       vir_loc(iv)%site=js
     endif
  enddo 
 endif

 nv=iv

!*********************************DEBUG Print Virtual List Local****************************
!*      write(*,*) 'Print Virtual List Local (IX, IY, Site)',ix,iy,is
!*      if(nv.gt.0) then
!*      do iv=1,nv
!*         write(*,*) vir_loc(iv)%x, vir_loc(iv)%y, vir_loc(iv)%site
!*      enddo
!*      endif
!********************************************************************************************
 RETURN
 END SUBROUTINE search_virtual_local

 SUBROUTINE virtual_status_change(ix,iy,is,vir_init,nv_init,vir_cre,n_cre,vir_anh,n_anh)
  USE calc_ads_en, only : constraint_satisfied
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ix,iy,is,nv_init
  INTEGER :: jx,jy,js                                             ! other site counter
  INTEGER :: ii,jj,kk                                             ! directions
  INTEGER :: ll,mm,nn                                             ! Some other counters
  INTEGER :: id                                                   ! index to count number of dependents.
  INTEGER :: i_cre                                                 ! index to count creation of virtual sites
  INTEGER :: n_cre, n_anh                                         ! number of created and annihilated sites.
  TYPE(DEPENDENT),INTENT(IN)::vir_init(2*max_dependent)           ! initial virtual sites
  TYPE(DEPENDENT) :: vir_cre(2*max_dependent)                     ! created virtual sites
  TYPE(DEPENDENT) :: vir_anh(2*max_dependent)                     ! annihilated virtual sites
  LOGICAL :: match_found

  i_cre=0
  n_cre=0
  n_anh=nv_init

  if(nv_init.eq.0) then
    call search_virtual_local(ix,iy,is,n_cre,vir_cre)             ! if previously nothing was in the list changes will be the only creation
!**         write(*,*) 'NV_INIT= 0'
    RETURN
  endif

!**         write(*,*) 'NV_INIT=',nv_init

  do ll=1,nv_init
  vir_anh(ll)=vir_init(ll)
  enddo 

  if(.not.constraint_satisfied(ix,iy,is,confg(ix,iy,is))) then
    match_found=.FALSE.
    do ll=1,nv_init 
      if(.not.match_found) then     
        if((ix.eq.vir_anh(ll)%x).and.(iy.eq.vir_anh(ll)%y).and.(is.eq.vir_anh(ll)%site)) then
           match_found=.TRUE.
!***                write(*,*) 'Match Found for',ix,iy,is

       endif

     else
        vir_anh(ll-1)=vir_anh(ll)
     endif   
   enddo

   if(match_found) then
     n_anh=nv_init-1                                         !If match found dimension of the list is reduced by 1                   
!           write(*,*) n_anh
    else
      n_anh=nv_init
      i_cre=i_cre+1
      vir_cre(i_cre)%x=ix
      vir_cre(i_cre)%y=iy
      vir_cre(i_cre)%site=is
    endif

  endif

  if(no_of_dependent(is).gt.0) then

    do id = 1,no_of_dependent(is)
     ii = dependent_list(is,id)%x
     jj = dependent_list(is,id)%y
     js = dependent_list(is,id)%site

     call site_finder(ix,iy,ii,jj,jx,jy)
     if(.not.constraint_satisfied(jx,jy,js,confg(jx,jy,js))) then
       match_found=.FALSE.
       if(n_anh.gt.0) then
        do ll=1,n_anh
         if(.not.match_found) then
           if((jx.eq.vir_anh(ll)%x).and.(jy.eq.vir_anh(ll)%y).and.(js.eq.vir_anh(ll)%site)) then
             match_found=.TRUE.
           endif
       
         else
          vir_anh(ll-1)=vir_anh(ll)
         endif
        enddo
       endif

       if(match_found) then
         n_anh=n_anh-1
       else
         i_cre=i_cre+1
         vir_cre(i_cre)%x=jx
         vir_cre(i_cre)%y=jy
         vir_cre(i_cre)%site=js
       endif
     endif
    enddo

  endif
  n_cre=i_cre

!************************DEBUG************************************
!***       if(n_anh.gt.0) then
!***         write(*,*) "Eureka! Annihilation of chaos!!!"
!***       endif
!***********************DEBUG*************************************

 RETURN
 END SUBROUTINE virtual_status_change

 SUBROUTINE update_virtual_list(n_cre,vir_cre,n_anh,vir_anh)
 IMPLICIT NONE
 INTEGER :: n_cre,n_anh       !number of virtual sites created or annihilated due to the move
 TYPE(DEPENDENT) :: vir_cre(2*max_dependent)
 TYPE(DEPENDENT) :: vir_anh(2*max_dependent)
 TYPE(DEPENDENT) :: temp_virtual_sites(0:nvirtual)

 INTEGER :: ii,jj,kk                                      !general counters
 INTEGER :: imatch,idump                                 !counting number of match found
 LOGICAL :: match                                         ! true if match found
 LOGICAL :: nvirtual_anh                                  ! number of virtual_sites after annihilation

 idump=0
 nvirtual_anh=nvirtual

 if((n_cre.eq.0).and.(n_anh.eq.0)) then
   RETURN
 endif

 if((n_anh.gt.0).and.(nvirtual.gt.0)) then 

  imatch=0
!Removing the annihilated virtual sites from the list
  do ii=1,nvirtual
   match=.FALSE.
   if(imatch.lt.n_anh) then
    do jj=1,n_anh     
     if((virtual_sites(ii)%x.eq.vir_anh(jj)%x).and.(virtual_sites(ii)%y.eq.vir_anh(jj)%y).and. &
         & (virtual_sites(ii)%site.eq.vir_anh(jj)%site)) then
        imatch=imatch+1
        match=.TRUE.
        exit
     endif
    enddo
   endif

   if(match) then
     cycle
   else
     idump=idump+1
     temp_virtual_sites(idump)= virtual_sites(ii)
   endif
  enddo
  
  nvirtual_anh=nvirtual_anh-n_anh                         ! no of virtual sites updated after annihilation

  do ii=1,nvirtual_anh                                    !virual_sites are updated after annihilation.
    virtual_sites(ii)=temp_virtual_sites(ii) 
  enddo  
 endif

 nvirtual=nvirtual-n_anh+n_cre

 if(n_cre.gt.0) then
   do ii=1,n_cre
      virtual_sites(nvirtual_anh+ii)%x=vir_cre(ii)%x
      virtual_sites(nvirtual_anh+ii)%y=vir_cre(ii)%y
      virtual_sites(nvirtual_anh+ii)%site=vir_cre(ii)%site
   enddo
 endif

 RETURN
 END SUBROUTINE update_virtual_list
 

 SUBROUTINE virtual_kawasaki(n_cre1,vir_cre1,n_cre2,vir_cre2,n_anh1,vir_anh1,n_anh2,vir_anh2,n_cre,vir_cre,n_anh,vir_anh)

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n_cre1,n_cre2,n_anh1,n_anh2
 TYPE(DEPENDENT),DIMENSION(2*max_dependent),INTENT(IN) :: vir_cre1,vir_cre2,vir_anh1,vir_anh2
 INTEGER :: n_cre, n_anh
 TYPE(DEPENDENT),DIMENSION(2*max_dependent) :: vir_cre, vir_anh
 
 if((n_cre1.eq.0).and.(n_cre2.eq.0)) then 
   n_cre=0
   goto   2120
 else if ((n_cre1.eq.1).and.(n_cre2.eq.0)) then
   vir_cre(1)=vir_cre1(1)
   n_cre=1
 else if ((n_cre1.eq.0).and.(n_cre2.eq.1)) then
   vir_cre(1)=vir_cre2(1)
   n_cre=1
 else
   call remove_redundancy(n_cre1,vir_cre1,n_cre2,vir_cre2,n_cre,vir_cre) 
 endif

2120  continue

 if((n_anh1.eq.0).and.(n_anh2.eq.0)) then
   n_anh=0
   goto   2121
 else if ((n_anh1.eq.1).and.(n_anh2.eq.0)) then
   vir_anh(1)=vir_anh1(1)
   n_anh=1
 else if ((n_anh1.eq.0).and.(n_anh2.eq.1)) then
   vir_anh(1)=vir_anh2(1)
   n_anh=1
 else
   call remove_redundancy(n_anh1,vir_anh1,n_anh2,vir_anh2,n_anh,vir_anh)
 endif

2121  continue
 RETURN
 END SUBROUTINE virtual_kawasaki

 SUBROUTINE remove_redundancy(na,a,nb,b,nc,c)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: na,nb
 INTEGER :: nc
 TYPE(DEPENDENT),DIMENSION(2*max_dependent),INTENT(IN) :: a,b
 TYPE(DEPENDENT) :: c(2*max_dependent)
     
 INTEGER :: ii,jj,idump
 LOGICAL :: match

 idump=0
 do ii=1,na
  match=.FALSE.
  do jj=1,nb  
     if((a(ii)%x.eq.b(jj)%x).and.(a(ii)%y.eq.b(jj)%y).and.(a(ii)%site.eq.b(jj)%site)) then
        match=.TRUE.
        exit
     endif
   enddo

   if(match) then
     cycle
   else
     idump=idump+1
     c(idump)%x=a(ii)%x
     c(idump)%y=a(ii)%y
     c(idump)%site=a(ii)%site
   endif
 enddo

 do ii=1,nb
    idump=idump+1
    c(idump)%x=b(ii)%x
    c(idump)%y=b(ii)%y
    c(idump)%site=b(ii)%site
 enddo

 nc=idump
 RETURN
 END SUBROUTINE remove_redundancy   


 FUNCTION kronecker_delta(n,m)
 INTEGER kronecker_delta
 INTEGER, INTENT(IN) :: n,m
 IF(n.eq.m) then
 kronecker_delta=1
 ELSE
 kronecker_delta=0
 ENDIF
 RETURN
 END FUNCTION

 END MODULE gcmc_module

      
