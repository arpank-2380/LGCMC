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


 SUBROUTINE read_simulation_control
 IMPLICIT NONE
!     Input related variables
 CHARACTER*100  :: buffer, label,label2,input
 INTEGER        :: posstart,posend  ! To search the key word label
 INTEGER        :: poseq,poscomment  ! position of "=" sign needed to search input  and position of comment flag
 INTEGER        :: ios=0, line=0
 INTEGER        :: count_component=0, count_site_type=0  !Counting the number of component and site type entered    
 INTEGER        :: count_temperature=0
 INTEGER        :: component_info_read_count=0
 INTEGER        :: corx,cory                             !Correlation function cut-off

 CHARACTER*80   :: component_temp_info, site_type_temp      !For reading temporary component name and site_type  
 INTEGER        :: compvar_temp, site_type_var_temp         !For reading temporary component and site_type variable
 INTEGER        :: comp_info_pos_strt, comp_info_pos_end    ! For reading component information: name
 CHARACTER*20   :: comp_info_label, comp_info_input(3)         ! Temporary variable to read comp_info
 INTEGER        :: comp_info_ios, comp_info_counter            ! takes care of input output status while reading comp_info from scratch file, counter variable
 INTEGER        :: comma_pos1, comma_pos2                   ! Commap ositions in component information

 DOUBLE PRECISION :: comp_info_data                            ! Temporary variable to store actual data while reading comp info
 LOGICAL        :: disp_start_read=.FALSE.               !whether displacement_start read or not
 LOGICAL        :: maxdisp_read=.FALSE.
 LOGICAL        :: corx_read=.FALSE., cory_read=.FALSE.
 LOGICAL        :: cons_en_cut_read = .FALSE.
 LOGICAL        :: file_exists                           !whether the SIMULATION_CONTROL file exists
 LOGICAL        :: comp_info_required = .TRUE.           !Whether component info (Tc, Pc, acc_fac) is required?        
 LOGICAL        :: comp_info_read_tmp(3),comp_info_read  ! Checking during reading each component info         

!     General counter related variables
 INTEGER        :: ii, icomp, jcomp

!ios is negative if an end of record condition is encountered or if an endfile condition is detected.
!ios is positive if an error was detected
!ios is zero otherwise

!all input related files are started with number 1 (101-199)
!all output related files are started with number 2 (201-299)
!all scratch files are started with number 3 (301-399) 


 INQUIRE(file='SIMULATION_CONTROL',exist=file_exists)

 if(file_exists) then
   open(101,file='SIMULATION_CONTROL',status='old')
 else
   stop_flag=.TRUE.
   write(*,'(A)') "!!!ERROR!!! SIMULATION_CONTROL file not found! Program will abort!"
   write(*,'(A)') "Please provide SIMULATION_CONTROL file to continue with the simulations"
   RETURN
 endif

 open(301, status='scratch')    !scratch file for temporary storage of component name
 open(302,status='scratch')    !scratch file for temporary storage of sity_type name

 Do While (ios.eq.0)
    
    label2=' '
    
    read(101,'(A)',iostat=ios) buffer
    if (ios.eq.0) then
       line=line+1
       poscomment=scan(buffer,'#!')
       if(poscomment.ne.0) then
         buffer=buffer(1:poscomment-1)
       endif
       
       posstart=verify(buffer,'   ')

       if (posstart.eq.0) then
          cycle
       endif
       
       posend = scan(buffer(posstart:),'- = ')+posstart-2
       poseq = scan(buffer,'=',.TRUE.)
       
!****            write(*,*) poscomment,posstart,posend,poseq
     
            
!c              print *, "WARNING!!! A blank line in SIMULATION_CONTROL!! at line no",line
            

       if (poseq.eq.0) then
         write(*,'(A)') "ERROR IN INPUT!!!! SIMULATION CONTROL. '=' sign NOT FOUND  at line no "//trim(int2str(line))

       else
          
!****            write(*,*)  poscomment,posstart,posend,poseq

         label= buffer(posstart:posend)
         label2=buffer(posend+2:poseq-1)
         input=buffer(poseq+1:)
        
         select case (label)
         case ("Temperature")
           read(input,*,iostat=ios) t
           count_temperature=count_temperature+1

         case ("Eqbm_step") 
           read(input,*,iostat=ios) neqbmstep

         case("Simu_step") 
           read(input,*,iostat=ios) nsimulstep

         case("EOS")
           read(input,*,iostat=ios) eos
 
         case("N_X_site")
           read(input,*,iostat=ios) nxsite

         case("N_Y_site")
           read(input,*,iostat=ios) nysite

         case("N_site_type")
           read(input,*,iostat=ios) nsite_type

         case("Restart")
           read(input,*,iostat=ios) restart 

         case("Corr_Func")
           read(input,*,iostat=ios) corr_func

         case("Corr_X")
           read(input,*,iostat=ios) corx
           if(ios.eq.0) then
             corx_read=.TRUE.
           endif

         case("Corr_Y")
           read(input,*,iostat=ios) cory
           if(ios.eq.0) then
             cory_read=.TRUE.
           endif
        
         case("Max_Cons")
           read(input,*,iostat=ios) max_cons

         case("N_component")
           read(input,*,iostat=ios) ncomp

         case("Lat_Int_On")
           read(input,*,iostat=ios) lat_int

         case("Gas_Int_On")
           read(input,*,iostat=ios) gas_int

         case("X_Cut")
           read(input,*,iostat=ios) xcut

         case("Y_Cut")
           read(input,*,iostat=ios) ycut
       
         case("Disp_Move_On")
           read(input,*,iostat=ios) displacement

         case("Disp_Move_Start")
           read(input,*,iostat=ios) disp_start
           if(ios.eq.0) then
             disp_start_read=.TRUE.
           endif

         case("Max_Disp")
           read(input,*,iostat=ios) maxdisp
           if(ios.eq.0) then
             maxdisp_read=.TRUE.
           endif       

         case("X_Disp")
           read(input,*,iostat=ios) x_disp
           
         case("Y_Disp")
           read(input,*,iostat=ios) y_disp 

         case ("Max_Trial_Disp")
            read(input,*,iostat=ios) maxtrial 

         case("Dump_Step")
           read(input,*,iostat=ios) tdump

         case("I_Seed")
           read(input,*,iostat=ios) iseed

         case("Cons_En_Cut")
            read(input,*,iostat=ios) h_cons_min
            if(ios.eq.0) then
              cons_en_cut_read=.TRUE.
            endif

         case ("Max_Dependent")
            read(input,*,iostat=ios) max_dependent

         case("Component")
           if(verify(label2,' ').eq.0) then
           write(*,'(A)') "!!!ERROR IN SIMULATION_CONTROL!!! Please input the COMPONENT INDEX &
       &(INTEGER) at line "//trim(int2str(line))
           stop_flag=.TRUE.
           endif

           read(label2,*,iostat=ios) compvar_temp
           write(301,'(I6)', advance = 'no') compvar_temp
           write(301,'(A)') input 
           count_component=count_component+1

         case("Site")
           if(verify(label2,' ').eq.0) then
           write(*,'(A)')"!!!ERROR IN SIMULATION_CONTROL!!! Please input the SITE INDEX (INTEGER) at line "//trim(int2str(line))
           stop_flag=.TRUE.
           endif
           
           read(label2,*,iostat=ios) site_type_var_temp
           read(input,*,iostat=ios)  site_type_temp
           write(302,*) site_type_var_temp,site_type_temp
!****                print *, site_type_var_temp,site_type_temp
           count_site_type=count_site_type+1

         
         case default
           write(*,'(A)') "Skipping invalid label of SIMULATION_CONTROL at line "//trim(int2str(line))
    
         end select  
       endif     
    endif
  Enddo

  ntotalstep=neqbmstep+nsimulstep
  nsite=nxsite*nysite    



! Printing error messages
 
     if (count_temperature.eq.0) then
        write(*,'(A)'),"!!!ERROR IN SIMULATION_CONTROL!!!  Please input a simulation Temperature"
     stop_flag=.TRUE.
     end if

     if (count_component.eq.0) then
      write(*,'(A)')"ERROR IN SIMULATION CONTROL!!! No component info"
      stop_flag=.TRUE.

     else if(count_component.ne.ncomp) then
      write(*,'(A)') "!!!ERROR IN SIMULATION_CONTROL!!! N_Component and number of supplied Gas Component mismatch"
      stop_flag=.TRUE.
     endif   
     
     if (count_site_type.eq.0) then
       write(*,'(A)')"!!!ERROR IN SIMULATION CONTROL!!! No site_type info"
       stop_flag=.TRUE.

     else if(count_site_type.ne.nsite_type) then
        write(*,'(A)') "!!!ERROR IN SIMULATION_CONTROL!!! Number of site information provided mismatch with N_site_type"
        stop_flag=.TRUE.
     endif

!Allocating arrays
 ALLOCATE(component(1:ncomp),site_type(1:nsite_type))
 ALLOCATE(comp_information(1:ncomp)) 


! Checking whether EOS is implemented

 select case (eos)
 case ("Ideal")
   eos_implemented=.TRUE.
   comp_info_required = .FALSE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Ideal gas law."

 case ("VdW")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Vander Waal's EOS."
 case ("RK")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Redlich-Kwong EOS."
   write(*,'(A)') "See Chem. Rev., 44, 233-244 (1949) for more details about RK EOS."

 case ("SRK")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Soave-Redlich-Kwong EOS."
   write(*,'(A)') "See Chem. Engg. Sci, 27, 1197-1203 (1972) for more details about SRK EOS."

 case ("PR")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Peng-Robinson EOS."
   write(*,'(A)') "See Industrial and Engineering Chemistry: Fundamentals, 15, 59-64 (1976) for more details about PR EOS."

 case ("PR78")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)')"Gas chemical potential will be calculated according to Peng-Robinson EOS modified in 1978."
   write(*,'(A)') "See Fluid Phase Equilibria, 447, 39-71 (2017); Table-1, Peng & Robinson (1978), for PR78 EOS"

 case ("PR80")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)')"Gas chemical potential will be calculated according to Peng-Robinson EOS modified in 1980"
   write(*,'(A)') "See Fluid Phase Equilibria, 447, 39-71 (2017); Table-1, Peng & Robinson (1980), for PR80 EOS"

 case ("PRG")
   eos_implemented=.TRUE.
   comp_info_required = .TRUE.
   write(*,'(A)') "Gas chemical potential will be calculated according to Peng-Robinson-Gasem EOS."
   write(*,'(A)') "See Fluid Phase Equilibria, 181, 113-125 (2001) for more details about PRG EOS."

 case default
   write(*,'(A)') "!!!WARNING!!! Entered EOS = "//trim(eos)//" is not recognized/implemented."
   write(*,'(A)') "You have to provide fugacities in GAS file."
   comp_info_required = .FALSE.

 end select



!Reading from scratch files

! Reading component informations

 REWIND (301) 

 ios=0
 Do while (ios.eq.0)
    read(301,'(I6)',iostat=ios,advance='no') compvar_temp
    read(301,'(A)',iostat=ios) component_temp_info

    if(ios.eq.0) then
      if ((compvar_temp.le.0).or.(compvar_temp.gt.ncomp)) then
        stop_flag=.TRUE.
        write(*,'(A)') "!!!ERROR IN SIMULATION CONTROL!!! Component-index do not lie between 1 and "//trim(int2str(ncomp))
        cycle
      endif

      comp_info_pos_strt = verify(component_temp_info,"    ")  ! extracting component name: start index from input 
      component_temp_info = component_temp_info(comp_info_pos_strt:)
      comp_info_pos_end = scan(component_temp_info," {")-1     ! extracting component name: end index from input 
      component(compvar_temp) = component_temp_info(1:comp_info_pos_end)  !storing component name

!****      write(8,'(A)') component_temp_info

         comp_info_pos_strt = index(component_temp_info,'{')
     
         if (comp_info_pos_strt.eq.0) then
            if (comp_info_required) then
               stop_flag=.TRUE.
               write(*,'(A)') "!!!ERROR IN SIMULATION_CONTROL!!! More information for component: "//trim(component(compvar_temp))//" is needed. See next line."
               write(*,'(A)') "Useage: Component-"//trim(int2str(compvar_temp))//" = "//trim(component(compvar_temp))//&
                          & " {Pc <FLOAT(in atm)>, Tc <FLOAT(in K)>, acc_fac <FLOAT> }"
            endif    
         else if (comp_info_required) then   
              comp_info_pos_end = index(component_temp_info,'}')
              component_temp_info = component_temp_info(comp_info_pos_strt+1:comp_info_pos_end-1)
              comma_pos1=scan(component_temp_info,',')
              comma_pos2=scan(component_temp_info,',',.TRUE.)
              comp_info_input(1)=component_temp_info(1:comma_pos1-1)
              comp_info_input(2)=component_temp_info(comma_pos1+1:comma_pos2-1)
              comp_info_input(3)=component_temp_info(comma_pos2+1:)
              
              comp_info_read_tmp(:)=(/.FALSE.,.FALSE.,.FALSE./)              
              do comp_info_counter=1,3
                 read(comp_info_input(comp_info_counter), *, iostat=comp_info_ios) comp_info_label, comp_info_data

!***** debug                 write(9,*) comp_info_label, comp_info_data              

                 if (comp_info_ios.eq.0) then
                    select case (comp_info_label) 
                           case("Pc")
                               comp_information(compvar_temp)%pc=comp_info_data
                               comp_info_read_tmp(1)=.TRUE.
                           case("Tc")
                               comp_information(compvar_temp)%tc=comp_info_data
                               comp_info_read_tmp(2)=.TRUE.
                           case("acc_fac")
                               comp_information(compvar_temp)%acc_fac=comp_info_data
                               comp_info_read_tmp(3)=.TRUE.
                    end select
                 endif
              enddo
            
              comp_info_read=comp_info_read_tmp(1).and.comp_info_read_tmp(2).and.comp_info_read_tmp(3) 
              if(.not.comp_info_read) then
                write(*,'(A)') "!!!ERROR in SIMULATION_CONTROL!!! Problem reading component info for "//trim(component(compvar_temp))//". See next line for suggestion."
                write(*,'(A)') "Useage: Component-"//trim(int2str(compvar_temp))//" = "//trim(component(compvar_temp))//&
                          & " {Pc <FLOAT(in atm)>, Tc <FLOAT(in K)>, acc_fac <FLOAT> }"

                stop_flag=.TRUE.
              endif 
          endif 
!****      write(8,'(A)') component_temp_info

    endif
 enddo

!**** Printing gas info for debug

! write(8,*)                                    "Component      Pc     Tc       acc_fac"
! Do comp_info_counter=1,ncomp
! write(8,'(A10,2X,G12.6,2x,G12.6,2x,G12.6)') component(comp_info_counter),comp_information(comp_info_counter)%pc,&
!                                           & comp_information(comp_info_counter)%tc, comp_information(comp_info_counter)%acc_fac 
! Enddo
!**** Debug blosk ends

 REWIND (302)

 ios=0
 Do while  (ios.eq.0)
    read(302,*,iostat=ios) site_type_var_temp,site_type_temp
    if(ios.eq.0) then
       if ((site_type_var_temp.le.0).or.(site_type_var_temp.gt.nsite_type)) then
          stop_flag=.TRUE.
          write(*,'(A)') "!!!ERROR IN SIMULATION CONTROL!!! Site-index do not lie between 1 and"//trim(int2str(nsite_type))
          cycle
       endif

       site_type(site_type_var_temp)=site_type_temp
   endif 
 enddo

close(101)
!Closing the scratch files
     close(301)
     close(302)

! if(eos.ne.'Ideal') then
!   write(*,'(A)') "Non-ideal equations of state (EOS) encountered!"
!   write(*,'(A)') "For each component, you have to provide fugacity coefficients according to "//trim(eos)//" EOS in gas.dat file." 
! endif

 if(disp_start_read) then
    if(disp_start.ge.neqbmstep) then
      write(*,'(A)') "ERROR!!! displacement move starting in production phase"
      write(*,'(A)') " Changing the start of displacement move at default value"
      disp_start=neqbmstep/2
    endif
 else
      disp_start=neqbmstep/2
 endif
 
 if(.not.maxdisp_read) then
   maxdisp= 2*nsite*nsite_type
 endif

  disp_box=(2*x_disp+1)*(2*y_disp+1)*nsite_type-1
 

   if(corx_read) then
     if(corx.gt.(nxsite/2)) then
       write(*,'(A)') "WARNING!!!Corr_X greater than 1/2 box length!!! Reseting Corr_X value to 1/2 box length"
       corx=nxsite/2
     endif
   else 
      corx=min(5,nxsite/2)
   endif

   if(cory_read) then
     if(cory.gt.(nysite/2)) then
       write(*,'(A)') "WARNING!!!Corr_Y greater the 1/2 box length!!!Reseting Corr_Y to 1/2 box length"
       cory=nysite/2
     endif
   else
     cory=min(5,nysite/2)
   endif

   if(.not.cons_en_cut_read) then
     h_cons_min=40*8.31*t/1000.0d0                 !scaled with temp-erature
   endif

   corx1=-1*corx
   corx2=corx
   cory1=-1*cory
   cory2=cory  


! Checking and Reading GAS_INT and GAS_INT_B files
 ALLOCATE(int_param_a(1:ncomp,1:ncomp),int_param_b(1:ncomp,1:ncomp))

!Initializing the interaction parameter to be 0
 do icomp=1,ncomp
    do jcomp=1,ncomp
       int_param_a(icomp,jcomp)=0.0d0
       int_param_b(icomp,jcomp)=0.0d0
    enddo
 enddo 


 if(eos_implemented.and.(eos.ne."Ideal").and.(ncomp.gt.1)) then
    write(*,'(A)', advance='no') "EOS mixing rule: "
    call read_gas_int (.TRUE.,'GAS_INT_A',106,gas_int_a_exists)
!    call read_gas_int (.FALSE.,'GAS_INT_B',107,gas_int_b_exists)
    if(gas_int_a_exists) then
      write(*,'(A)') "a[i,j] = (1-k[i,j])*sqrt(a[i]*a[j]); k[i,i] = 0 always but k[i,j] = 0 unless provided in GAS_INT_A."
    else
      write(*,'(A)') "a[i,j] = sqrt(a[i]*a[j])"
    endif

    call read_gas_int (.FALSE.,'GAS_INT_B',107,gas_int_b_exists)    
    write(*,'(17X)',advance='no')
    if(gas_int_b_exists) then
      write(*,'(A)') "b[i,j] = (1-l[i,j])*0.5*(b[i]+b[j]); l[i,i] = 0 always but l[i,j] = 0 unless provided in GAS_INT_B."
    else
      write(*,'(A)') "b[i,j] = 0.5*(b[i]+b[j])"
    endif
 endif

! If user mistakenly make int_param_a(i,i) !=0 or int_param_b(i,i) != 0 we
! correct for it
  do icomp = 1, ncomp
     int_param_a(icomp,icomp) = 0.0d0
     int_param_b(icomp,icomp) = 0.0d0
  enddo



!Printing
!*           Print *, "Temperature =",t
!*           Print *, "Eqbm_step =",neqbmstep
!*           Print *, "Simu_step =",nsimulstep
!*           Print *, "Total_step =",ntotalstep
!*           Print *, "N_X_site =", nxsite
!*           Print *, "N_Y_Site =", nysite
!*           Print *, "N_Site =", nsite
!*           Print *, "N_site_type", nsite_type
!*           Print *, "N_comp =", ncomp
!*           Print *, "EOS = ", eos
!*           Print *, "Restart = ", restart
!*           do ii=1,ncomp
!*           print *, ii,component(ii)
!*           enddo           
!*           do ii=1,nsite_type
!*           print *, ii, site_type(ii)
!*           enddo
!      STOP
!      END
!do icomp=1,ncomp-1
!   do jcomp=icomp+1,ncomp
!      write(*,*) trim(component(icomp)),' ', trim(component(jcomp)), int_param_a(icomp,jcomp),int_param_b(icomp,jcomp)
!   enddo
!enddo

 RETURN
 END SUBROUTINE read_simulation_control
