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


 SUBROUTINE read_lateral_interactions
 IMPLICIT NONE
 INTEGER, PARAMETER :: min_input=7      

 DOUBLE PRECISION, PARAMETER :: epsj=0.00010d0
 DOUBLE PRECISION, PARAMETER :: sym_no=2.0d0
!     Variable needed to read input file
 CHARACTER*200  :: full_input_line,buffer
 CHARACTER*40   :: label(min_input+1), input(min_input+1)
 INTEGER        :: posstart,posend ! To search the key word label
 INTEGER        :: poscomment ! position of comment flag 
 INTEGER        :: ios=0, line=0, ilabel=0, i_filled_line=0 ! i_filled_line = index of non-empty line which is read by the program
 INTEGER        :: nlabel,ninput ! total number of label(key-word) and input specified
 INTEGER        :: count_lat_int_read=0 ! counting number of adsorption energy read by the program
      
!     counter variable 
 INTEGER        :: icomp1,icomp2,isite_type1,isite_type2 !component and site type index of two sites
 INTEGER        :: ix,iy   !vector shift between two sites
 INTEGER        :: ii,jj,kk
 
 
 DOUBLE PRECISION :: lat_int_temp !for temporary storage of lateral interaction
 LOGICAL        :: component_mismatch_flag =.FALSE. ! TRUE when a mismatch of component is found
 LOGICAL        :: site_type_mismatch_flag =.FALSE. ! TRUE when a mismatch of site type is found
!***      LOGICAL,ALLOCATABLE :: read_lat_int(:,:,:,:,:,:) ! Whether the array element has been read
 LOGICAL        :: read_line_status
!     logical variables to keep track reading
 LOGICAL        :: read_comp1,read_comp2,read_site1,read_site2
 LOGICAL        :: read_ix,read_iy,read_lat
 LOGICAL        :: asym_temp,asym_found
 LOGICAL        :: file_exists                        !To check whether LAT_INT file exist

 INTEGER :: xmin, ymin, xmax, ymax

!     Variable to check inconsistency
 LOGICAL, ALLOCATABLE :: print_error(:,:,:,:,:,:) !Needed while checking similarities of different lateral interactions and whether lattice assymmetry needed
 LOGICAL              :: print_error_here
 LOGICAL, ALLOCATABLE :: chk_done(:,:,:,:,:,:)    !Whether data inconsistency checking has been done or not!!!

 LOGICAL              :: expression(6)        !needed to check inconsistency of data when isite_type1=isite_type2

 xmin=-1*xcut
 xmax=xcut
 ymin=-1*ycut
 ymax=ycut
 
      
!     Allocating lateral interaction energy array
 ALLOCATE(j_raw(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))

 ALLOCATE(j(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))

 ALLOCATE(jb(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))
 
!     Allocating book-keeping array
 ALLOCATE(read_lat_int(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))

!     Allocating lattice-assymmetry array
 ALLOCATE(asym(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))

!     Allocating printing error_message needed array
 ALLOCATE(print_error(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))

!     Allocating checking done array
 ALLOCATE(chk_done(1:ncomp,1:nsite_type,-1*xcut:xcut,-1*ycut:ycut,1:ncomp,1:nsite_type))
      
!     Initializing allocated arrays
      
 j(:,:,:,:,:,:)=0.0d0
 
 read_lat_int(:,:,:,:,:,:)=.FALSE.

 asym(:,:,:,:,:,:)=.FALSE.

 print_error(:,:,:,:,:,:)=.TRUE.
 
 chk_done(:,:,:,:,:,:)=.FALSE.
!     Opening files
    
 INQUIRE(file='LAT_INT',exist=file_exists)

 if(file_exists) then    
   Open(103,file='LAT_INT',status='old')
 else
   stop_flag=.TRUE.
   write(*,'(A)') "!!!ERROR!!! LAT_INT file not found! Program will abort!"
   write(*,'(A)') "Either set Lat_Int_On = F or provide LAT_INT file to continue with the simulations"
   RETURN
 endif
 
 do while (ios.eq.0)

!     Initializing logical arrays to keep track reading  as FALSE
    read_comp1=.FALSE.
    read_comp2=.FALSE.
    read_site1=.FALSE.
    read_site2=.FALSE.
    read_ix=.FALSE.
    read_iy=.FALSE.
    read_lat=.FALSE.
    asym_temp=.FALSE.

    read(103,'(A)',iostat=ios) full_input_line
    
    if(ios.eq.0) then
       line=line+1
       buffer=full_input_line
       poscomment=scan(buffer,'!#')

!****      print *,line, full_input_line
      
       if (poscomment.ne.0) then
          buffer=buffer(1:poscomment-1)
       endif
       posstart=verify(buffer,'   ')
       
       if (posstart.eq.0) then
!*     print *, "WARNING!!! Skipping blank line in ADS_FREE_EN in",line
          cycle
       endif
            
!     Reading labels (key-words)
            
       i_filled_line=i_filled_line+1
       
       if (i_filled_line.eq.1) then
          
          do while((posstart.ne.0).and.(ilabel.lt.(min_input+1)))
             ilabel=ilabel+1
             posend=scan(buffer(posstart:),'   ')+posstart-2
             label(ilabel)=buffer(posstart:posend)
             buffer=buffer(posend+1:)
             posstart=verify(buffer,'   ')
          enddo
          
          nlabel=ilabel
               
               
               
!****  write(*,*) nlabel,(label(ilabel),ilabel=1,nlabel)
               
               
!     Reading inputs 
               
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
          
          if (ninput.lt.min_input) then
             stop_flag=.TRUE.
             write(*,'(A)')' less number of input than required at line '//trim(int2str(line))//' in LAT_INT'
             write(*,'(A)') 'At least 6 input (Comp1,Site1,IX,IY,Comp2,Site2,Lat_Int) needed' 
          endif
          
          do ilabel=1,ninput
             select case (label(ilabel))   
             case('Component1')
!****  write(*,*) 'Component found',input(ilabel)
!****  write(*,*) integer_test(input(ilabel))                     
                if(integer_test(input(ilabel))) then
                   read(input(ilabel),*,iostat=ios) ii
                   if((ii.gt.0).and.(ii.le.ncomp)) then
                      icomp1=ii
                      component_mismatch_flag=.FALSE.
                      read_comp1=.TRUE.
!****  write(*,*) ii
                   else
                      component_mismatch_flag=.TRUE.
                   endif
                   
                else
                   do ii=1,ncomp
                      if(component(ii).eq.input(ilabel)) then
                         icomp1=ii
                         component_mismatch_flag=.FALSE.
                         read_comp1=.TRUE.
!****  write(*,*) ii
                         exit
                      else
                         component_mismatch_flag=.TRUE.
                      endif
                   enddo
                endif

                if(component_mismatch_flag) then
                  write(*,'(A)') "ERROR in LAT_INT!!! Component1 mismatch found at line "//trim(int2str(line))
                  stop_flag=.TRUE.
                endif
                
             case('Component2')
!****  write(*,*) 'Component found',input(ilabel)
!****  write(*,*) integer_test(input(ilabel))                     
                if(integer_test(input(ilabel))) then
                   read(input(ilabel),*,iostat=ios) ii
                   if((ii.gt.0).and.(ii.le.ncomp)) then
                      icomp2=ii
                      component_mismatch_flag=.FALSE.
                      read_comp2=.TRUE.
!****  write(*,*) ii
                   else
                      component_mismatch_flag=.TRUE.
                   endif
                   
                else
                   do ii=1,ncomp
                      if(component(ii).eq.input(ilabel)) then
                         icomp2=ii
                         component_mismatch_flag=.FALSE.
                         read_comp2=.TRUE.
!****  write(*,*) ii
                         exit
                      else
                         component_mismatch_flag=.TRUE.
                      endif
                   enddo
                endif

                if(component_mismatch_flag) then
                   write(*,'(A)') "ERROR in LAT_INT!!! Component2 mismatch found at line "//trim(int2str(line))
                  stop_flag=.TRUE.
                endif

                
             case('Site1')
!****  write(*,*) 'Site found: ', input(ilabel)
                if(integer_test(input(ilabel))) then
!****  write(*,*) integer_test(input(ilabel)) 
                   read(input(ilabel),*,iostat=ios) ii
                   if((ii.gt.0).and.(ii.le.nsite_type)) then
                      isite_type1=ii
                      site_type_mismatch_flag=.FALSE.
                      read_site1=.TRUE.
!****  write(*,*) isite_type
                   else
                      site_type_mismatch_flag=.TRUE.
                   endif
                   
                else
                   do ii=1,nsite_type
                      if(site_type(ii).eq.input(ilabel)) then
                         isite_type1=ii
                         site_type_mismatch_flag=.FALSE.
                         read_site1=.TRUE.
!****  write(*,*) isite_type
                         exit
                      else
                         site_type_mismatch_flag=.TRUE.
                      endif
                   enddo
                endif

                if(site_type_mismatch_flag) then
                  write(*,'(A)') "ERROR in LAT_INT!!! Site1 mismatch found at line "//trim(int2str(line))
                  stop_flag=.TRUE.
                endif

                
             case('Site2')
!****  write(*,*) 'Site found: ', input(ilabel)
                if(integer_test(input(ilabel))) then
!****  write(*,*) integer_test(input(ilabel)) 
                   read(input(ilabel),*,iostat=ios) ii
                   if((ii.gt.0).and.(ii.le.nsite_type)) then
                      isite_type2=ii
                      site_type_mismatch_flag=.FALSE.
                      read_site2=.TRUE.
!****  write(*,*) isite_type
                   else
                      site_type_mismatch_flag=.TRUE.
                   endif
                   
                else
                   do ii=1,nsite_type
                      if(site_type(ii).eq.input(ilabel)) then
                         isite_type2=ii
                         site_type_mismatch_flag=.FALSE.
                         read_site2=.TRUE.
!****  write(*,*) isite_type
                         exit
                      else
                         site_type_mismatch_flag=.TRUE.
                      endif
                   enddo
                endif
                
                if(site_type_mismatch_flag) then
                  write(*,'(A)') "ERROR in LAT_INT!!! Site2 mismatch found at line "//trim(int2str(line))
                  stop_flag=.TRUE.
                endif

                
             case('IX')
                read(input(ilabel),*,iostat=ios) ii
                if(ios.ne.0) then
                   write(*,'(A)')'Error in LAT_INT!!! IX should be an integer between (' //trim(int2str(xmin))//':'&
                &  //trim(int2str(xmax))//') at line '//trim(int2str(line))
                   stop_flag=.TRUE.
                else if((ii.lt.(-1.0*xcut).or.(ii.gt.xcut))) then
                   stop_flag=.TRUE.
                   write(*,'(A)') 'Error in LAT_INT!!! IX out of bound at line'// trim(int2str(line))//'in LAT_INT'
                   Print*,'IX should be within (-XCUT:XCUT)'
                else
                   ix=ii
                   read_ix=.TRUE.
                endif
                
             case('IY')
                read(input(ilabel),*,iostat=ios) ii
                if(ios.ne.0) then
                   write(*,'(A)')'Error in LAT_INT!!! IY should be an integer between ('//trim(int2str(ymin))//':'& 
                &  //trim(int2str(ymax))//') at line '//trim(int2str(line)) 
                   stop_flag=.TRUE.
                else if ((ii.lt.(-1.0*ycut).or.(ii.gt.ycut))) then
                   stop_flag=.TRUE.
                   write(*,'(A)'),'Error in LAT_INT!!! IY out of bound at line'//trim(int2str(line))//' in LAT_INT'
                   Print*,'IY should be within (-YCUT:YCUT)'
                else
                   iy=ii
                   read_iy=.TRUE.
                endif

             case('Lat_Int')
                   read(input(ilabel),*,iostat=ios) lat_int_temp
                   if(ios.eq.0) then
                      read_lat=.TRUE.
                   endif
             
             case('Lattice_Asym')
                 read(input(ilabel),*,iostat=ios) asym_temp 

                   
             case default
                   Print*,'Error in LAT_INT!!! Keyword not recognised'
                   stop_flag=.TRUE.
                   exit
             end select
          enddo

                
          read_line_status=read_comp1.and.read_comp2.and.read_site1.and.read_site2.and.read_ix.and.read_iy.and.read_lat
          
          if(read_line_status) then
             if((ix.eq.0).and.(iy.eq.0).and. &       !If both the chosen sites are same 
            &     (isite_type1.eq.isite_type2)) then !Lateral interaction has to be 0 = initialized value
                write(*,'(A)') 'Warning!!!Chosen sites are same!!! There could not be any lateral interaction' 
                ! This is removal of self interaction error if used do a mistake
                write(*,'(A)'), 'Neglecting entered data and self interaction will set to be  0'
                read_lat_int(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.FALSE.
                
                go to 1004
                
             else if(read_lat_int(icomp1,isite_type1,ix,iy,icomp2,isite_type2)) then
                write(*,'(A)')'Warning!!! Data repeated in LAT_INT at line'//trim(int2str(line))
                Print*,'Data with higher line index will be considered'
             endif
             
             j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=lat_int_temp

             
             read_lat_int(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
       
             asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=asym_temp
                  
 1004        continue

          else
             write(*,'(A)') "ERROR in LAT_INT!!! Incomplete Line at line"//trim(int2str(line))              
             stop_flag=.TRUE.   
          endif
       endif
    endif
 enddo 
      
!     File reading is complete now
!     Processing the lateral_interaction matrix j
!     h(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=h(icomp2,isite_type2,-ix,-iy,icomp1,icomp2)
      
 Do icomp1=1,ncomp
    Do icomp2=1,ncomp
       Do isite_type1=1,nsite_type
          Do isite_type2=1,nsite_type
             Do ix=xmin,xmax
                Do iy=ymin,ymax
                   
                   if((ix.eq.0).and.(iy.eq.0).and.(isite_type1.eq.isite_type2)) then
                      go to 1005
                   else if(read_lat_int(icomp1,isite_type1,ix,iy,icomp2,isite_type2)) then
                      if(.not.read_lat_int(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)) then 

!Lateral interaction is unique between 2 specified sites. between them any one could be working site and the other one could be neighbor site
!So exchange of label (icomp1,isite_type1) with (icomp2,isite_type2) and (ix,iy) with (-ix,-iy) is defining same lateral interaction.

                         j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)
                         
                         read_lat_int(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE. 
!This flag is set true now as this data is now stored in array.  
                      endif                                       
                                              
                   endif
                           
1005               continue
                   enddo
                enddo     
             enddo    
          enddo
     enddo
 enddo      
              

! Processing lateral_interaction matrix j assumption: If exchange of component between 2 site do not change lateral interaction
! But if the interaction is provided in the LAT_INT file it will take care of that.
 Do icomp1=1,ncomp
    Do icomp2=1,ncomp
       Do isite_type1=1,nsite_type
          Do isite_type2=1,nsite_type
             Do ix=xmin,xmax
                Do iy=ymin,ymax

                   if((ix.eq.0).and.(iy.eq.0).and.(isite_type1.eq.isite_type2)) then
                        cycle
                   endif
                   
                   if(.not.read_lat_int(icomp2,isite_type1,ix,iy,icomp1,isite_type2)) then

                         j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2) 

!!!assumption: If exchange of component between 2 site do not change lateral interaction 
!!! If lateral interaction provided in the LAT_INT file then it will skip this assumption.

                     read_lat_int(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=.TRUE.
! This flag is set true now as this data is now stored in array. 
                   endif
               
                enddo
              enddo
           enddo
        enddo
     enddo
  enddo

              
! Initializing jb array

 Do icomp1=1,ncomp
    Do icomp2=1,ncomp
       Do isite_type1=1,nsite_type
          Do isite_type2=1,nsite_type
             Do ix=xmin,xmax
                Do iy=ymin,ymax
                
                 jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)= j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
           
                Enddo
             Enddo
          Enddo
       Enddo
    Enddo
 Enddo


! Checking whether everything is consistent or not
       
 Do icomp1=1,ncomp
    Do icomp2=1,ncomp
       Do isite_type1=1,nsite_type
          Do isite_type2=1,nsite_type
             Do ix=xmin,xmax
                Do iy=ymin,ymax

!****             asym_found=asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)

                     print_error_here=print_error(icomp1,isite_type1,ix,iy,icomp2,isite_type2)

                      if (chk_done(icomp1,isite_type1,ix,iy,icomp2,isite_type2)) then
                        go to 1006

                      else if((ix.eq.0).and.(iy.eq.0).and.(isite_type1.eq.isite_type2)) then
                         go to 1006

                      else if (icomp1.eq.icomp2) then
                         if(Abs(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)-j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)).gt.epsj) then 

                           if((ix.eq.0).and.(iy.eq.0)) then
                              stop_flag=.TRUE.
                              if(print_error_here) then            !print_error_if1
                              write(*,*) "Error in LAT_INT!!!"
                              write(*,*) "Assymmetry of interaction is not allowed when both IX and IY are 0 together"
                              write(*,*) "This violates the periodicity of the system"
                              write(*,*) 

                              write(*,*) "See more info below!!!"
                              endif                                !print_error_endif1
                           else
                              if(print_error_here) then            !print_error_if2
                              write(*,*) "WARNING!!!"
                              endif                                !print_error_endif2
                           endif
         
                           

                           if(print_error_here) then               !print_error_if3
                           write(*,'(A)') "Assymmetry of interacti on has been detected between following neighbours of component "// trim(component(icomp1)) 
                           write(*,'(A)') "Site1= "//trim(site_type(isite_type1))//" Site2= "//trim(site_type(isite_type2))//" IX= " &
                        &   //trim(int2str(ix))//" IY= "//trim(int2str(iy)) 
                           write(*,'(A)')"Lateral interaction energy for this pair of site/component "//&
                        &  trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2 )))//'  '//trim(dble2str(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)))

                           write(*,*) 
                           endif                                   !print_error_endif_3
                           
                           asym_found=asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                        
                           
                           if((ix.ne.0).or.(iy.ne.0)) then       
                              if(isite_type1.eq.isite_type2) then
                                 if(.not.asym_found) then
                                    if (print_error_here) then     !print_error_if4
                                    write(*,*) "WARNING!!!"
                                    write(*,'(A)') "No Lattice Assymmetry info was found for Site "//trim(site_type(isite_type1))// " and component "//trim(component(icomp1)) 
                                    write(*,'(A)') "But for same component and same site type assymmetry of interaction is only possible by Lattice Assymmetry"

                                    write(*,'(A)') "WARNING!!! "//"Imposing lattice Assymmetry for component "//trim(component(icomp1))//" and site "//trim(site_type(isite_type1))
                                    write(*,*) 
                                    endif                          !print_error_endif4

                                 else
                                   write(*,'(A)') "WARNING!!!"
                                   write(*,'(A)') "Lattice Assymmetry was found in LAT_INT for component "//trim(component(icomp1))// &
                                 & " between sites "//trim(site_type(isite_type1))//" and "//trim(site_type(isite_type2))
                                   write(*,'(A)') "Imposing Lattice Assymmetry for component "//trim(component(icomp1))//" between sites "//& 
                                 & trim(site_type(isite_type1))//" and "//trim(site_type(isite_type2))
                                   write(*,*)

                                endif
                             
                                    asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.

                                    asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.
                                    
                                
                                   jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)

                                   jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)

                              else if (asym_found) then
                                 if(print_error_here) then         !print_error_if5
                                 write(*,'(A)') "WARNING!!!"
                                 write(*,'(A)') "Lattice Assymmetry was found in LAT_INT for component "//trim(component(icomp1))//" between sites "//& 
                               &  trim(site_type(isite_type1))//" and "//trim(site_type(isite_type2))
                                 write(*,'(A)') "Imposing Lattice Assymmetry for component "//trim(component(icomp1))//" between sites "// & 
                               &  trim(site_type(isite_type1))//" and "//trim(site_type(isite_type2)) 
       
                                 write(*,*) 
                                 endif                             !print_error_endif5

                                    asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
                                    asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.

                                    jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)= j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)

                                    jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)

                              endif

                           if(print_error_here) then               !print_error_if6
                           write(*,'(A)') "Left hand and right hand coupling will be different along direction IX= "//trim(int2str(ix))//" IY= "//trim(int2str(iy))
                           write(*,'(100A)') ('-',kk=1,100)
                           write(*,*)
                           endif                                   !print_error_endif6
                           endif
                         
                           print_error(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.FALSE.

                         else
                            if((asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1))) then
                               if(print_error_here) then          !print_error_if6.0.1
                                write(*,'(A)')"WARNING about ERROR!!!"
                                write(*,'(A)') "Lattice Assymmetry flag is ON!!! for Comp1=Comp2="//trim(component(icomp1))// &
                              &   " Site1="// trim(site_type(isite_type1))//" and Site2="//trim(site_type(isite_type2))// & 
                              &   " along IX="//trim(int2str(ix))//" and IY="// trim(int2str(iy))
                                write(*,'(A)') "But no assymmetry of interaction has been detected along this direction"
                                write(*,'(A)') "If you want to invoke the Lattice Assymmetry please supply another Lat_Int parameter"
                                write(*,'(A)') "Removing ERROR!!!"
                                write(*,'(A)') "Revoking Lattice Assymmetry"
                                write(*,'(100A)') ('-',kk=1,100)
                                write(*,*)
                                endif                             !print_error_if6.0.1
                                asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2) = .FALSE.
                                asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .FALSE.
                  
                            endif

                         endif

                           chk_done(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.
                           chk_done(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.

                              
                                    
                      else if (isite_type1.eq.isite_type2) then
!Checking inconsistency of data:
                        expression(1)=Abs(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)-j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)).gt.epsj   
                        ! if it returns .TRUE. that means they are not same

                        expression(2)=Abs(j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)-j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)).gt.epsj

                        expression(3)=Abs(j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)-j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)).gt.epsj

                        expression(4)=Abs(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)-j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)).gt.epsj

                                             
                        expression(5)=Abs(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)-j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)).gt.epsj
                        ! Diagonal expression-1. needed to print the error message

                        expression(6)=Abs(j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)-j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)).gt.epsj
                        !Diagonal expression-2


                           asym_found=asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2).or.  &
                      &               asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)



                        jj=0
                        do ii=1,4
                           if(.not.expression(ii)) then
                             jj=jj+1
                           endif
                        enddo 
                    
                        if((jj.ne.2).and.(jj.ne.4)) then
                           stop_flag=.TRUE.
                           if(print_error_here) then               !print_error_if6.1
                              write(*,'(A)') "ERROR in LAT_INT!!!" 
                              write(*,'(A)') "More than 2 Lat_Int parameter found for site "//trim(site_type(isite_type1))//" for Comp1="// & 
                            & trim(component(icomp1))//" and Comp2="//trim(component(icomp2))//"along direction IX="//trim(int2str(ix))// & 
                            & " and IY="// trim(int2str(iy))
                              write(*,'(A)')"Due to symmetry and periodicity of the system, it can't have more that 2 numbers"

                              write(*,'(A)') "Incompatible Lat_Int's are as follows:"
                              write(*,'(A)') "Comp1="//trim(component(icomp1))//" Comp2="//trim(component(icomp2))//" IX="//trim(int2str(ix))// & 
                            &  " IY="//trim(int2str(iy))//" Lat_Int="//trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)))

                              if(expression(1)) then           ! Are this 2 number different then print
                                 write(*,'(A)') "Comp1="//trim(component(icomp2))//" Comp2="//trim(component(icomp1))//" IX="//trim(int2str(ix))// & 
                            &    " IY="//trim(int2str(iy))//" Lat_Int="//trim(dble2str(j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)))
                              endif

                              if(expression(2).and.expression(5)) then
                                 write(*,'(A)') "Comp1="//trim(component(icomp1))//" Comp2="//trim(component(icomp2))//" IX="// trim(int2str(-1*ix))// & 
                                &  " IY="//trim(int2str(-1*iy))//" Lat_Int="//trim(dble2str(j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type2)))
                              endif

                              if(expression(3).and.expression(4).and.expression(6)) then
                                 write(*,'(A)') "Comp1="//trim(component(icomp2))//" Comp2="//trim(component(icomp1))//" IX="//trim(int2str(-1*ix))// & 
                          &     " IY="//trim(int2str(-1*iy))//" Lat_Int="//trim(dble2str(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)))
                              endif

                              write(*,'(A)') "Please correct the Error to run the simulation!!!"
                              write(*,'(100A)')('-',kk=1,100)
                              write(*,*) 
                           endif                                   !print_error_endif6.1

                        else if(jj.eq.4) then
                          if(asym_found) then
                            if(print_error_here) then              !print_error_if6.1.1
                            write(*,'(A)')"WARNING about ERROR!!!"
                            write(*,'(A)') "Lattice Assymmetry tag is ON!!! for Comp1=Comp2="//trim(component(icomp1))//" Site1="// & 
                        &    trim(site_type(isite_type1))//" and Site2="//trim(site_type(isite_type2))//" along IX="//trim(int2str(ix))//" and IY="//trim(int2str(iy))
                            write(*,'(A)') "But no Assymmetry of interaction has been detected"
                            write(*,'(A)') "If you want to introduce Lattice Assymmetry supply another Lat_Int parameter!"
                            write(*,'(A)') "Removing ERROR!!"
                            write(*,'(A)') "Lattice Assymmetry tag will be neglected during the simulation"
                            write(*,'(A)') "Revoking Lattice Assymmetry!!!"
                            write(*,'(100A)') ('-',kk=1,100)
                            write(*,*)
                            endif                                  !print_error_endif6.1.1
                          endif
              
                          
                             asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2) = .FALSE.
                             asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2) = .FALSE.
                             asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1) = .FALSE.
                             asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .FALSE.


                             chk_done(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
                             chk_done(icomp2,isite_type1,ix,iy,icomp1,isite_type2)= .TRUE.
                             chk_done(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)= .TRUE.
                             chk_done(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .TRUE.

                             go to 1006                           
                        endif

!****                                asym_found= asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2).or. &
!****                              &             asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)


                        if((.not.expression(1)).and.expression(2).and.(.not.expression(3)).and.expression(4)) then

                           if(print_error_here) then               !print_error_if6.2
                           write(*,*) "WARNING!!!"
                           write(*,'(A)')"Assymmetry of interaction has been detected between following neighbours of site "//trim(site_type(isite_type1))
                           write(*,'(A)') " Comp1= "//trim(component(icomp1))//" Comp2= "//trim( component(icomp2))//" IX= "//trim(int2str(ix))//" IY= "//trim(int2str(iy))

                           write(*,'(A)')"Lateral interaction energy for this pair of site/component "//trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)))// &
                       &   '  '//trim(dble2str(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)))
                           write(*,*)

                           write(*,'(A)') "Parameters specified for site="//trim(site_type(isite_type1))//" along the direction IX="//trim(int2str(ix))// & 
                       &   "and IY="//trim(int2str(iy))//" for Comp1="//trim(component(icomp1))//" and Comp2="//trim(component(icomp2))
                           write(*,'(A)') "can be valid only if there are !!!Lattice Assymmetry!!! along this direction for this set of components  for that site"
                           if(.not.asym_found) then
                             write(*,'(A)') "WARNING!!!"
                             write(*,'(A)') " Lattice Assymmetry not found in LAT_INT. But it will be imposed"
                           endif

                           write(*,'(A)') " Imposing Lattice Assymmetry"
                           write(*,'(100A)') ('-',kk=1,100) 
                           write(*,*)
                           endif                                   !print_error_endif6.2

                           asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2) = .TRUE.
                           asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2) =.TRUE.
                           asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1) = .TRUE.
                           asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .TRUE.


                           jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                           jb(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type2)
                           jb(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)=j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)
                           jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)

                        else if(expression(1).and.(.not.expression(2)).and.expression(3).and.(.not.expression(4))) then

                           if(print_error_here) then               !print_error_if6.3
                           write(*,*) "WARNING!!!"
                           write(*,'(A)')"Assymmetry of interaction has been detected between following neighbours of site "//trim(site_type(isite_type1))
                           write(*,'(A)') " Comp1= "//trim(component(icomp1))//" Comp2= "//trim( component(icomp2))//" IX= "//trim(int2str(ix))//" IY= "//trim(int2str(iy))

                           write(*,'(A)') "Lateral interaction energy for this pair of site/component "//trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)))// & 
                        &  '  '//trim(dble2str(j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)))
                           write(*,*)
                           write(*,'(A)') "Parameters specified for site="//trim(site_type(isite_type1))//" along the direction IX=" & 
                        &  //trim(int2str(ix))//"and IY="//trim(int2str(iy))//" for Comp1="//trim(component(icomp1))//" and Comp2="//trim(component(icomp2))
                           write(*,'(A)') "can be valid only if there are !!!NO Lattice Assymmetry!!! along this direction for this set of components  for that site"

                           if(asym_found) then
                             write(*,'(A)') "WARNING and ERROR!!!"
                             write(*,'(A)') "Lattice Assymmetry has been detected for this pair of component/site along this direction"
                             write(*,'(A)') "You cannot have lattice Assymmetry according to the defined lateral interactions!!!!"
!****                        stop_flag=.TRUE.
                             write(*,'(A)') "REMOVING ERROR!!!"
                             write(*,'(A)') "Neglecting the Lattice Assymmetry tag!!!"
                             write(*,'(A)') "Simulation will proceed without any Lattice Assymmetry for this pair of components along this direction"
                             write(*,'(A)') "Revoking Lattice Assymmetry"
                           endif
                           write(*,'(100A)')('-',kk=1,100)
                           write(*,*)
                           endif                                  !print_error_endif6.4
                           
                           asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.FALSE.
                           asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2) =.FALSE.
                           asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1) = .FALSE.
                           asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .FALSE.



                        else if(expression(5).neqv.expression(6)) then

                          if(print_error_here) then                !print_error_if7
                          write(*,*) "WARNING!!!"
                          write(*,'(A)') "Assymmetry of interaction has been detected between following neighbours of site "//trim(site_type(isite_type1))
                          write(*,'(A)') " Comp1= "//trim(component(icomp1))//" Comp2= "//trim( component(icomp2))//" IX= "//trim(int2str(ix))//" IY= "//trim(int2str(iy))

                          if(expression(5)) then
                          write(*,'(A)') "Lateral interaction energy for this pair of site/component "//trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)))// & 
                     &                   '  '//trim(dble2str(j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)))

                          else
                          write(*,'(A)')"Lateral interaction energy for this pair of site/component "//trim(dble2str(j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)))// &
                     &      " "//trim(dble2str(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)))
                          endif
                          write(*,*)
                          endif                                    !print_error_endif7
                         

                          if(asym_found) then
                             if(print_error_here) then             !print_error_if8
                             write(*,'(A)') "WARNING!!!"
                             write(*,'(A)') "Lattice assymmetry was found in LAT_INT for site "//trim(site_type(isite_type1))//" between components "// &
                     &          trim(component(icomp1))//" and "// trim(component(icomp2))
                             write(*,'(A)') "Imposing Lattice Assymmetry"  

                             write(*,*)
                             endif                                 !print_error_endif8
                             asym(icomp1,isite_type1,ix,iy,icomp2, isite_type2)=.TRUE.
                             asym(icomp2,isite_type1,ix,iy,icomp1,isite_type2) =.TRUE.
                             asym(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1) = .TRUE.
                             asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .TRUE.

                             if(expression(5)) then

                             jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)
                             jb(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)
                             jb(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)
                             jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)

                             else if (expression(2)) then

                             jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)
                             jb(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)
                             jb(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)=j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)
                             jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp2,isite_type1,ix,iy,icomp1,isite_type2)

                             else  

                             jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                             jb(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                             jb(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                             jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)
                             endif

                          endif
                          
                          if(print_error_here) then                !print_error_if9
                          write(*,'(A)') "Left hand and right hand coupling will be different along direction IX= "//trim(int2str(ix))//" IY= "//trim(int2str(iy))
                          write(*,'(100A)') ('-',kk=1,100)
                          write(*,*)
                          endif                                    !print_error_endif9
                          print_error(icomp2,isite_type1,ix,iy,icomp1,isite_type2)=.FALSE.
                          print_error(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)=.FALSE.
                          print_error(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type2)=.FALSE.
                        endif


                        chk_done(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
                        chk_done(icomp2,isite_type1,ix,iy,icomp1,isite_type2)= .TRUE.
                        chk_done(icomp1,isite_type2,-1*ix,-1*iy,icomp2,isite_type1)= .TRUE.
                        chk_done(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1) = .TRUE.

                      else
                         if(Abs(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)-j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)).gt.epsj)  then
              
                            if(print_error_here) then              !print_error_if10
                            write(*,'(A)') "WARNING!!!"     
                            write(*,'(A)') "Assymmetry of interaction has been detected between Sites "//trim(site_type(isite_type1))//" and "//trim(site_type(isite_type2))//" and components "// &
                     &        trim(component(icomp1))//" and "//trim(component(icomp2))//" along direction IX="//trim(int2str(ix))//" IY="//trim(int2str(iy))

                            write(*,'(A)')"Lateral interaction energy for this pair of site/component "//trim(dble2str(j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)))//'  '// &
                     &       trim(dble2str(j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)))

                            endif                                  !print_error_endif10
                         
                            asym_found=asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2).or.asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)

                            if(.not.asym_found) then
                                    if(print_error_here) then      !print_error_if11
                                    write(*,*) "WARNING!!!"
                                    write(*,'(A)') "No Lattice Assymmetry info was found for Sites "//trim(site_type (isite_type1))// " and "//trim(site_type(isite_type2))// & 
                                 &    " and components "//trim(component(icomp1))//" and "//trim(component(icomp2))//" along direction IX="//trim(int2str(ix))//" IY="//trim(int2str(iy)) 
                                    write(*,'(A)') "But Site1="//trim(site_type(isite_type1))//" Comp1="//trim(component(icomp1))//" IX="//trim(int2str(ix))//" IY="// & 
                                 &    trim(int2str(iy))//" Site2="//trim(site_type(isite_type2))// "Comp2="//trim(component(icomp2))
                                    write(*,'(A)') "and"
                                    write(*,'(A)') "Site1="//trim(site_type(isite_type2))//" Comp1="//trim(component(icomp2))//" IX="//trim(int2str(-1*ix))//" IY="// & 
                                 &    trim(int2str(-1*iy))//" Site2="//trim(site_type(isite_type1))//"Comp2="//trim(component(icomp1)) 
                                    write(*,'(A)') "represent same interaction"
                                    write(*,'(A)') "They can only be different if there is Lattice Assymmetry"
                                    write(*,*)  
                                    write(*,'(A)') "Imposing Lattice assymmetry along this direction for above mentioned sites and coxmponents"
                                    write(*,'(100A)')('-',kk=1,100)
                                    write(*,*)
     
!****                                         asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
!****                                         asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.
                                    
                                    endif                          !print_error_endif11
                            else
                                    if(print_error_here) then      !print_error_if12
                                    write(*,*) "Lattice Assymmetry was found for sites "//trim(site_type(isite_type1))// " and "//trim(site_type(isite_type2))// & 
                                 &    " and components "//trim(component(icomp1))//" and "//trim(component(icomp2))//" along direction IX="//trim(int2str(ix))// & 
                                 &    " IY="//trim(int2str(iy))
                                    write(*,'(A)') "Imposing Lattice assymmetry along this direction for above mentioned sites and components"
                                    write(*,'(100A)')('-',kk=1,100)
                                    endif                          !print_error_endif12

                           endif 

                           asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=.TRUE.
                           asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.

                           jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=j(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)
                           jb(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)


                         else     
                             if(asym_found) then
                                if(print_error_here) then         !print_error_if13
                                write(*,'(A)') "WARNING about ERROR!!!"
                                write(*,'(A)') "Lattice Assymmetry tag is ON!!! for"
                                write(*,'(A)') "Comp1="//trim(component(icomp1))//" Site1="//trim(site_type(isite_type1))//" Comp2="//trim(component(icomp2))// &
                             &   " Site2="//trim(site_type(isite_type2))//" along IX="//trim(int2str(ix))//" IY="//trim(int2str(iy))
                                write(*,'(A)')"But no assymmetry of interaction has been detected"
                                write(*,'(A)')"If you want to have Lattice Assymmetry invoke another Lat_Int parameter"
                                write(*,'(A)') "Removing Error!!!!"
                                write(*,'(A)') "Revoking Lattice Assymmetry!!!"
                                write(*,'(100A)') ('-',kk=1,100)
                                write(*,*)
                                endif                              !print_error_endif14

                                asym(icomp1,isite_type1,ix,iy,icomp2,isite_type2)  =  .FALSE.
                                asym(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)  =.FALSE.

                             endif
                         endif

                         print_error(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.FALSE.
                         chk_done(icomp1,isite_type1,ix,iy,icomp2,isite_type2)= .TRUE.
                         chk_done(icomp2,isite_type2,-1*ix,-1*iy,icomp1,isite_type1)=.TRUE.

                           endif
1006                    continue

                 Enddo
              Enddo 
           Enddo
        Enddo
     Enddo
  Enddo
      
!Changing the j according to symmetry_number
 Do icomp1=1,ncomp
    Do icomp2=1,ncomp
       Do isite_type1=1,nsite_type
          Do isite_type2=1,nsite_type
             Do ix=xmin,xmax
                Do iy=ymin,ymax

                j_raw(icomp1,isite_type1,ix,iy,icomp2,isite_type2)= j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)



                lat_int_temp=j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)/sym_no

                j(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=lat_int_temp



                lat_int_temp=jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)/sym_no
                
                jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)=lat_int_temp


!****       write(*,*) icomp1,icomp2,isite_type1,isite_type2,ix,iy,j(icomp1,isite_type1,ix,iy,icomp2,isite_type2),jb(icomp1,isite_type1,ix,iy,icomp2,isite_type2)

                 Enddo
              Enddo
           Enddo
        Enddo
     Enddo
  Enddo

                     
 

!Deaallocating unnecessary arrays

 DEALLOCATE(chk_done,print_error)
    
 RETURN         
 END SUBROUTINE
