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



 SUBROUTINE read_ads_free_en

 IMPLICIT NONE

!Input related variables
 CHARACTER*100  :: full_input_line,buffer
 CHARACTER*40   :: label(5), input(5)
 INTEGER        :: posstart,posend  ! To search the key word label
 INTEGER        :: poscomment            ! position of comment flag 
 INTEGER        :: ios=0, line=0, ilabel=0, i_filled_line=0   ! i_filled_line = index of non-empty line which is read by the program
 INTEGER        :: nlabel,ninput                              ! total number of label(key-word) and input specified
 INTEGER        :: count_ads_en_read=0                        ! counting number of adsorption energy read by the program
 INTEGER        :: icomp=1,isite_type=1                       ! index of component and site_type
 INTEGER        :: ads_en_to_read                             ! number of adsorption energy data to read  
!Temporary storage variables
 DOUBLE PRECISION :: ads_free_en_temp                         ! for temporary storage of adsorption free energy 
 LOGICAL          :: constraint_present_temp                  ! constraint_present

 LOGICAL        :: constraint_flag_read=.FALSE.                ! Determine whether Constraints tag given in File ADS_FREE_EN
 LOGICAL        :: component_mismatch_flag =.FALSE.           ! TRUE when a mismatch of component is found
 LOGICAL        :: site_type_mismatch_flag =.FALSE.           ! TRUE when a mismatch of site type is found
 LOGICAL        :: read_icomp=.FALSE., read_isite_type=.FALSE.     
 LOGICAL        :: file_exists                                !TRUE when ADS_FREE_EN file exists
!Declaring counter variables
 INTEGER ii,jj     

!Declaration of allocatable arrays needed if "Site" or "Component" Key-word missing
 INTEGER,ALLOCATABLE :: read_count_comp(:)                   ! Needed if "Site" Keyword missing. keep track the icomp entered
 INTEGER,ALLOCATABLE :: read_count_site_type(:)              ! Needed if "Component" keyword missing. keep track the isite_type
 INTEGER,ALLOCATABLE :: read_count_comp_site(:,:)            ! Needed to remove redaundancy error

!Allocating adsorption energy and constraints_present
 ALLOCATE(h(1:ncomp,1:nsite_type))
 ALLOCATE(constraint_present(1:ncomp,1:nsite))

! ALLOCATING  optional read arrays needed for default case
 ALLOCATE(read_count_comp(1:ncomp))
 ALLOCATE(read_count_site_type(1:nsite_type))
 ALLOCATE(read_count_comp_site(1:ncomp,1:nsite_type))       ! to check redaundancy error 

!Initializing constraint_present arrays with false ! Default is no constraint
 constraint_present(:,:) = .FALSE.

!Initializing read_count_arrays
 read_count_comp(:)=0
 read_count_site_type(:)=0
 read_count_comp_site(:,:)=0


 ads_en_to_read=ncomp*nsite_type   !No of adsorption energies to read base on SIMULATION_CONTROL

!Opening files
  
 INQUIRE(file='ADS_FREE_EN',exist=file_exists)
 if (file_exists) then
    Open(102,file='ADS_FREE_EN',status='old')
 else
    stop_flag=.TRUE.
    write(*,'(A)') "!!!ERROR!!! ADS_FREE_EN file not found! Program will abort!"
    write(*,'(A)') "Please provide ADS_FREE_EN file to continue with the simulations."
    RETURN
 endif
 

 do while (ios.eq.0)
!Initializing constraint_present_temp as FALSE as default one is false
    constraint_present_temp=.FALSE.
    read(102,'(A)',iostat=ios) full_input_line

    if(ios.eq.0) then
      line=line+1
      buffer=full_input_line
      poscomment=scan(buffer,'!#')
      if (poscomment.ne.0) then
         buffer=buffer(1:poscomment-1)
      endif 
      posstart=verify(buffer,'   ')   

      if (posstart.eq.0) then
!*             print *, "WARNING!!! Skipping blank line in ADS_FREE_EN in"
!*     x  ,line
      cycle         
      endif


!        Reading labels (key-words)

      i_filled_line=i_filled_line+1

      if (i_filled_line.eq.1) then 

       do while((posstart.ne.0).and.(ilabel.lt.5))
         ilabel=ilabel+1
         posend=scan(buffer(posstart:),'   ')+posstart-2
         label(ilabel)=buffer(posstart:posend)
         buffer=buffer(posend+1:)
         posstart=verify(buffer,'   ') 
       enddo
    
       nlabel=ilabel 

!****            write(*,*) nlabel,(label(ilabel),ilabel=1,nlabel)


!       Reading inputs 

      else if (i_filled_line.gt.1) then
       ilabel=0
       do while ((posstart.ne.0).and.(ilabel.lt.nlabel))
         ilabel=ilabel+1
         posend=scan(buffer(posstart:),'   ')+posstart-2
         input(ilabel)=buffer(posstart:posend)
         buffer=buffer(posend+1:)
         posstart=verify(buffer,'   ')
       enddo
       ninput=ilabel

       if (ninput.lt.nlabel) then
       write(*,'(A)')'WARNING! less number of input in line '//trim(int2str(line))//' than label specified in ADS_FREE_EN'
       endif 
 
!****            write(*,*) (input(ilabel),ilabel=1,nlabel)
       do ilabel=1,ninput
         select case(label(ilabel))

         case('Component') 
           read_icomp=.TRUE.  
!****              write(*,*) 'Component found',input(ilabel)
!****              write(*,*) integer_test(input(ilabel))                     
           if(integer_test(input(ilabel))) then
             read(input(ilabel),*,iostat=ios) ii
             if((ii.gt.0).and.(ii.le.ncomp)) then
               icomp=ii
               component_mismatch_flag=.FALSE.                   
!****                    write(*,*) ii
             else
               component_mismatch_flag=.TRUE.
             endif 
     
           else
             do ii=1,ncomp 
               if(component(ii).eq.input(ilabel)) then
                icomp=ii
                component_mismatch_flag=.FALSE.
!****                     write(*,*) ii
                exit
               else 
                component_mismatch_flag=.TRUE.
               endif
             enddo 
           endif

         case('Site') 
           read_isite_type=.TRUE.
!****                write(*,*) 'Site found: ', input(ilabel)
           if(integer_test(input(ilabel))) then            
!****                  write(*,*) integer_test(input(ilabel)) 
             read(input(ilabel),*,iostat=ios) ii
             if((ii.gt.0).and.(ii.le.nsite_type)) then
               isite_type=ii
               site_type_mismatch_flag=.FALSE.
!****                    write(*,*) isite_type
             else
               site_type_mismatch_flag=.TRUE. 
             endif

           else
             do ii=1,nsite_type
               if(site_type(ii).eq.input(ilabel)) then
                 isite_type=ii
                 site_type_mismatch_flag=.FALSE.
!****                      write(*,*) isite_type
                 exit
               else
                site_type_mismatch_flag=.TRUE.
               endif
             enddo
           endif 

         case('Ads_Free_En') 
!****                write(*,*) "Ads_Free_En found:", input(ilabel)
           read(input(ilabel),*,iostat=ios) ads_free_en_temp

!****                write(*,*) ads_free_en_temp

         case('Constraints')
           constraint_flag_read=.TRUE.
           read(input(ilabel),*,iostat=ios) constraint_present_temp
!*                if(ios.lt.0) then
!*                  constraint_present_temp=.FALSE.
!*                  ios=0
!*                endif

         case default
           stop_flag=.TRUE.
           write(*,'(A)')"ERROR!!!Label not recognised in ADS_FREE_EN"
           RETURN
         end select
       enddo

! Cases if icomp and/or isite_type is not supplied in ADS_FREE_EN
       if(.not.(read_icomp.or.read_isite_type)) then
         if(i_filled_line.gt.2) then
           if(isite_type.lt. nsite_type) then
             isite_type=isite_type+1
           else if(icomp.lt.ncomp) then
             isite_type=1
             icomp=icomp+1
           else
             icomp=1
             isite_type=1
           endif
         endif

       else if (.not.(read_icomp)) then       ! it did not read icomp i.e it read isite_type
            read_count_site_type(isite_type)=read_count_site_type(isite_type)+1
            if(read_count_site_type(isite_type).le.ncomp) then
              icomp=read_count_site_type(isite_type)
            else
              icomp=1
              read_count_site_type(isite_type)=1
            endif

       else if (.not.(read_isite_type)) then   !it did not read isite_type i.e it read icomp
            read_count_comp(icomp)=read_count_comp(icomp)+1
            if(read_count_comp(icomp).le.nsite_type) then
              isite_type=read_count_comp(icomp)
            else
              isite_type=1
              read_count_comp(icomp)=1
            endif
       endif                                 

            read_count_comp_site(icomp,isite_type)= read_count_comp_site(icomp,isite_type)+1
       
     
!****           write(*,*) icomp, isite_type

        
! writing the input line in scratch file       

!****            write(*,*)  'Hello I reached line here'
!****            write(*,*) icomp,isite_type,ads_free_en_temp

! Feeding the data in adssorption free energy array and constraint_present array
       h(icomp,isite_type)=ads_free_en_temp

       if(constraint_flag_read) then
         constraint_present(icomp,isite_type)= constraint_present_temp
       endif

! Checking the error and taking necessary actions 
       if(site_type_mismatch_flag) then
         stop_flag=.TRUE.
         write(*,'(A)') "Site type mismatch found in ADS_FREE_EN at line "//trim(int2str(line))
       endif
       if  (component_mismatch_flag) then
         stop_flag=.TRUE.
         write(*,'(A)') "Component type mismatch found in ADS_FREE_EN at line "//trim(int2str(line))
       endif
         
       count_ads_en_read=count_ads_en_read+1 
!****            write(303,*) icomp,isite_type,ads_free_en_temp
       

     endif     

!*        else
!*          print *, 'File read complete'
!*          print *, 'line=',line,'read count',count_ads_en_read
   endif
 enddo         

!    Checking if there are any discrepency in the supplied data enen if there are sufficient line when both icomp and isite_type supplied

 do icomp=1,ncomp
   do isite_type=1,nsite_type 
    if(read_count_comp_site(icomp,isite_type).lt.1) then
      stop_flag=.TRUE.
      write(*,'(A)') 'No ADS_FREE_EN data for (component/site) '// trim(component(icomp)) // '/'//trim(site_type(isite_type))
      if(count_ads_en_read.ge.ads_en_to_read) then
        Print *, 'Please check for redundancy!'
      endif
    endif
   enddo
 enddo


!     Deallocation of the arrays needed to read default case when either 'Site' or 'Component' is missing
 DEALLOCATE(read_count_comp,read_count_site_type)
 DEALLOCATE(read_count_comp_site)
 
 if(count_ads_en_read.lt.ads_en_to_read) then
   stop_flag=.TRUE.
   Print *, "!!!ERROR!!! Insufficient Ads_Free_En data"
 else if (count_ads_en_read.gt.ads_en_to_read) then
   Print *, "!!!WARNING!!!Confused!!! More Ads_Free_En data than  required!!! &
   & In redaundant case data in higher line index of ADS_FREE_EN will be taken into account"
 go to 1001
 endif

1001  continue

 if(stop_flag) then
   Print *, "!!!ERROR in ADS_FREE_EN!!!"
   RETURN
 endif 
 

!!!      h(0,:)=0.0D0    !Adsorption energy of vacuum  just for convinience

!  Determining whether constraints file needed or not
 Do icomp=1,ncomp
  Do isite_type=1,nsite_type
    constraint_file_needed=constraint_file_needed.or.constraint_present(icomp,isite_type)
    if (constraint_file_needed) then 
       Print *, 'You have constraints on LG Hamiltonian'
       Print *, 'You need to provide a CONSTRAINTS file'
       goto 1002
    endif
  Enddo
 Enddo

1002  continue
      
 RETURN
 END SUBROUTINE read_ads_free_en  


