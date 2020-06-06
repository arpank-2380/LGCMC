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


 USE variable
 USE logo_module
 USE cubic_eos_module, ONLY: calculate_fugacity 
 USE read_input_files
 USE gcmc_module

 IMPLICIT NONE
 INTEGER icomp,ipoint,ix,iy,ievn_odd
 INTEGER icomp1,icomp2,isite_type1,isite_type2,i,is
 DOUBLE PRECISION x_adsorbed
 REAL t1,t2
 DOUBLE PRECISION, ALLOCATABLE :: avgocc(:,:),avgocc_evn_odd(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: dist_f(:,:,:,:,:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: corr_f(:,:,:,:,:,:,:)
!      DOUBLE PRECISION              :: d_t
!      DOUBLE PRECISION              :: c_t
     
 DOUBLE PRECISION :: x_sum
 DOUBLE PRECISION, ALLOCATABLE :: total_pop(:)

 INTEGER time
 CHARACTER*50 ctime
 CHARACTER*80 corr_func_path
 CHARACTER*4  corr_func_var
 CHARACTER*20 temp_comp
 INTEGER      cmnd_argu_count
 CHARACTER*4  chk_input
 CHARACTER*80,ALLOCATABLE :: fname(:)

  cmnd_argu_count=command_argument_count()  
  if(cmnd_argu_count.gt.0) then
    call get_command_argument(1,chk_input)
  endif

  call cpu_time(t1) 
  
  call credit_print(6)

  write(*,'(141A)') ('-',i=1,141)
  write(*,'(141A)')(' ',i=1,55), "Messages related to input files",(' ',i=1,55)
  write(*,'(141A)') ('-',i=1,141)

  call read_simulation_control

  if(.not.stop_flag) then
    call read_ads_free_en
  endif

  if((.not.stop_flag).and.constraint_file_needed) then
     call read_constraints
  endif

  if((.not.stop_flag).and.lat_int) then
    call read_lateral_interactions
  endif

  if(.not.stop_flag) then
    call read_gas_data
  endif

  write(*,'(141A)') (' ',i=1,50), ('*',i=1,41), (' ',i=1,50)
  if (stop_flag) then
    write(*,'(141A)') (' ',i=1,50),'* !!!ERROR in input files!!! Aborted!!! *' , (' ',i=1,50)
    write(*,'(141A)') (' ',i=1,50),('*',i=1,41),(' ',i=1,50)
    STOP  
  else 
    write(*,'(141A)') (' ',i=1,50), '*  Input Sucessful. Check the Warnings  *'
    write(*,'(141A)') (' ',i=1,50), ('*',i=1,41), (' ',i=1,50)
  endif
  

  if(chk_input .eq. "-chk") then
    go to 3000
  endif  

  if(.not.allocated(x)) then
    ALLOCATE(x(1:ncomp))
  endif

  if (.not.allocated(phi)) then
    ALLOCATE(phi(1:ncomp))
  endif

  if (.not.allocated(f)) then
    ALLOCATE(f(0:ncomp))
  endif

  if(.not.allocated(avgocc)) then
    ALLOCATE(avgocc(ncomp,nsite_type))
  endif

  if(.not.allocated(dist_f)) then
    ALLOCATE(dist_f(ncomp,nsite_type,corx1:corx2,cory1:cory2,ncomp,nsite_type,9))

  endif

  if(.not.allocated(corr_f)) then
    ALLOCATE(corr_f(ncomp,nsite_type,corx1:corx2,cory1:cory2,ncomp,nsite_type,9))
  endif

  if (.not.allocated(avgocc_evn_odd)) then
     ALLOCATE(avgocc_evn_odd(ncomp,nsite_type,4))
  endif

  if(.not.allocated(fname)) then
     ALLOCATE(fname(ncomp)) 
  endif

  if(.not.allocated(total_pop)) then
     ALLOCATE(total_pop(ncomp))
  endif

  do icomp=1,ncomp
     fname(icomp)='allP_'//trim(component(icomp))//'.dat'
  enddo

 !File opening

 if(ncomp.gt.1) then
 open(202,file="allP.dat",status="unknown") 
 call credit_print(202)   

 write(202,'(A)',advance='no') "#"
 write(202,17) "   P(atm)   ",('   y_'//trim(component(icomp))//'    ',icomp=1,ncomp-1), &
              & (' theta_total_'//trim(component(icomp)),icomp=1,ncomp)  

 write(202,'(A)') "#"
 endif

 do icomp=1,ncomp
    temp_comp=component(icomp)
    open(202+icomp,file=fname(icomp),status='unknown')
    call credit_print(202+icomp)
    write(202+icomp,'(A)',advance='no') "#"
    write(202+icomp,18) "   P(atm)   ","y_"//trim(temp_comp)//'    ',(' theta_'//trim(site_type(i))//&
                        & '_'//trim(temp_comp),i=1,nsite_type),' theta_total_'//trim(temp_comp)
    write(202+icomp,'(A)') "#"
 enddo

 if(corr_func) then
   call system ("mkdir corr_func")
 endif
 
 write(*,*)  
 write(*,*) "Starting job at  ",ctime(time())
 write(*,*) 

 do ipoint=1,npoint
   write(*,*) 
   write(*,*) 
   write(*,'(A)') "----------------------------------"
   write(*,'(A)') "| Starting a new GCMC simulation |"              
   write(*,'(A)') "----------------------------------"

   p=gas_p(ipoint)
   do icomp=1,ncomp
      x(icomp)=gas_comp(ipoint,icomp)
      if(.not.eos_implemented) then
        phi(icomp)=gas_fug(ipoint,icomp)
      endif
   enddo
   if(eos_implemented) then
     call calculate_fugacity(p,t,x,phi)
   endif


   write(*,'(A)') "|**** Gas Information ****|"

   write(*,'(1X,A7,G14.6)',advance='no') "P(atm)=",p
   do icomp=1,ncomp-1
      write(*,'(2X,A12,1X,G14.6)',advance='no') "  y_"//trim(component(icomp))//"= ",x(icomp)
   enddo
   write(*,'(2X,A12,1X,G14.6)') "  y_"//trim(component(ncomp))//"= ",x(ncomp)
   if(eos_implemented) then
     write(*,'(22A)', advance='no') (' ',i=1,22)   
     do icomp=1,ncomp-1
        write(*,'(2X,A12,1X,G14.6)',advance='no') "phi_"//trim(component(icomp))//"= ",phi(icomp)
     enddo
     write(*,'(2X,A12,1X,G14.6)') "phi_"//trim(component(ncomp))//"= ",phi(ncomp)
   endif
!   write(*,'(A)') "--------------------------------"
   write(*,'(A)') "|**** End of fugacity calculation ****|"
!   write(*,'(A)') "-------------------------------"



   write(*,*)
   write(*,*) "Starting GCMC at  ", ctime(time())

!   write(*,'(1X,A7,G14.6)',advance='no') "P(atm)=",p
!   do icomp=1,ncomp-1
!      write(*,'(2X,A9,1X,G14.6)',advance='no') "y_"//trim(component(icomp))//"= ",x(icomp)
!   enddo
!   write(*,'(2X,A9,1X,G14.6)') "y_"//trim(component(ncomp))//"=",x(ncomp)
!   if(eos_implemented) then
!     do icomp=1,ncomp-1
!        write(*,'(2X,A12,1X,G14.6)',advance='no') "phi_"//trim(component(icomp))//"= ",phi(icomp)
!     enddo
!     write(*,'(2X,A12,1X,G14.6)') "phi_"//trim(component(ncomp))//"= ",phi(ncomp)
!   endif

   call gcmc(avgocc,dist_f,corr_f,avgocc_evn_odd)
   
   do icomp=1,ncomp
     total_pop(icomp)=0.0d0
     do i=1,nsite_type
        total_pop(icomp)=total_pop(icomp)+avgocc(icomp,i)
     enddo             
     write(202+icomp,15) p,x(icomp),(avgocc(icomp,i),i=1,nsite_type),total_pop(icomp)
   enddo

   if (ncomp.gt.1) then
   write(202,16) p,(x(icomp),icomp=1,ncomp-1),(total_pop(icomp),icomp=1,ncomp)
   endif

    if(ncomp.eq.2) then
      if(nsite_type.eq.1) then
        x_adsorbed=avgocc(1,1)/(avgocc(1,1)+avgocc(2,1))
        write(202,10) p,x(1),avgocc(1,1),avgocc(2,1),x_adsorbed

      else if (nsite_type.eq.2) then
        x_adsorbed=(avgocc(1,1)+avgocc(1,2))/(avgocc(1,1)+avgocc(1,2)+avgocc(2,1)+avgocc(2,2))
        write(202,10) p,x(1),(avgocc(1,1)+avgocc(1,2)),avgocc(2,1)+avgocc(2,2),x_adsorbed
        write(3001,10) p,x(1),avgocc(1,1),avgocc(2,1),avgocc(1,1)/(avgocc(1,1)+avgocc(2,1))       ! What are these printing? Which file? Where is file with unit 3001 & 3002?
        write(3002,10) p,x(1),avgocc(1,2),avgocc(2,2),avgocc(1,2)/(avgocc(1,2)+avgocc(2,2))       ! What are these printing?
      endif

    else if(ncomp.eq.1) then
      if(nsite_type.eq.2) then
        write(202,11) p,avgocc(1,1),avgocc(1,2)
      else
        write(202,12) p,avgocc(1,1)
      endif
    endif

!Wriiting correlation functions
   if(corr_func) then
     write(corr_func_var,'(I4.4)') ipoint
     corr_func_path=trim('corr_func/')//trim("simulation_")//trim(corr_func_var)
     
     open(4000+ipoint,file=corr_func_path,status='unknown')
     call credit_print(4000+ipoint)

     write(4000+ipoint,*) "Pressure (atm) =",p
     do icomp=1,ncomp-1
      write(4000+ipoint,'(2X,A9,1X,G14.6)',advance='no') "y_"//trim(component(icomp))//"=",x(icomp)
     enddo
     write(4000+ipoint,'(2X,A9,1X,G14.6)') "y_"//trim(component(ncomp))//"=",x(ncomp)


     do icomp1=1,ncomp
        do icomp2=1,ncomp
           do isite_type1=1,nsite_type
              do isite_type2=1,nsite_type

                 write(4000+ipoint,'(100A)') ('*',i=1,100)
                 write(4000+ipoint,'(A)') "Correlation function between "//trim(component(icomp1))//"and"//trim(component(icomp2))
                 write(4000+ipoint,'(A)') "Site-1 is a "//trim(site_type(isite_type1))//" and Site-2 is a "//trim(site_type(isite_type2))
                 write(4000+ipoint,'(100A)') ('*',i=1,100)

                 do ievn_odd=1,9
                
                    if(ievn_odd.eq.1) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"X=EVEN and Y=EVEN"
                       write(4000+ipoint,'(100A)') ('x',i=1,50) 
                       
                    else if (ievn_odd.eq.2) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"X=EVEN and Y=ODD"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                    else if (ievn_odd.eq.3) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"X=ODD and Y=EVEN"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                    else if (ievn_odd.eq.4) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50) 
                       write(4000+ipoint,'(A)')"X=ODD and Y=ODD"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                    else if (ievn_odd.eq.5) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"X=EVEN"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                    else if (ievn_odd.eq.6) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"Y=EVEN"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                   else if (ievn_odd.eq.7) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"X=ODD"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                   else if (ievn_odd.eq.8) then
                       write(4000+ipoint,'(100A)') ('x',i=1,50)
                       write(4000+ipoint,'(A)')"Y=ODD"
                       write(4000+ipoint,'(100A)') ('x',i=1,50)

                   else 

                       write(4000+ipoint,'(100A)') ('-',i=1,100)
                       write(4000+ipoint,'(A)') "Total Correlation functions without Even odd discrimination"
                       write(4000+ipoint,'(100A)') ('-',i=1,100)

                   endif

 
                    write(4000+ipoint,'(A)') "  IX     IY   "//"      dist_f        "//"      corr_f       "

                    do ix=corx1,corx2
                       do iy=cory1,cory2                     

                          write(4000+ipoint,14) ix,iy,dist_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,ievn_odd), &
                                                & corr_f(icomp1,isite_type1,ix,iy,icomp2,isite_type2,ievn_odd)
                       enddo
                    enddo
                 enddo
              enddo   
           enddo
        enddo
     enddo      
   endif
 enddo 

10 Format(F9.5,2x,F9.5,2x,G15.8,2x,G15.8,2x,G15.8)
11 Format(F9.5,2x,G15.8,2x,G15.8)
12 Format(F9.5,2x,G15.8)
14 Format(I5,2x,I5,2X,G15.8,2X,G15.8)
15 Format(1X,G12.6,4x,G12.6,4X,<nsite_type+1>G20.8)
16 Format(1x,G12.6,4x,<ncomp-1>G15.6,4x,<ncomp>G20.8)
17 Format(A12,4x,<ncomp-1>A15,4x,<ncomp>A20)
18 Format(A12,4x,A12,4x,<nsite_type+1>A20)
 write(*,*) "GCMC finished at  ", ctime(time())


3000  continue
 call cpu_time(t2)
!****      write(*,*) (avgocc(icomp),icomp=1,2)
 write(*,*)
 write(*,*)
 write(*,*) "job finished in", t2-t1,"secs" 

 STOP
 END

