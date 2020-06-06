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

     
 SUBROUTINE read_constraints 
 IMPLICIT NONE


 CHARACTER*40, PARAMETER :: default_inp_fmt_label(1:5) =(/'Site','X','Y','Component','Present'/)

 CHARACTER*20, PARAMETER :: default_label(1:3)=(/'Component','Site','Cons_En'/)


 DOUBLE PRECISION :: cons_en_temp=0.0d0               !Temporary storage of constrained_free_energy
      
!Variables needed to read input file
 CHARACTER*200:: full_input_line,buffer
 CHARACTER*100:: comp_buff                           ! Needed to analyse component input
 CHARACTER*40 :: inp_fmt_label(1:5)
 INTEGER      :: ninp_fmt_label
 LOGICAL      :: default_inp_fmt=.TRUE.               ! Determine whether the input_format is default
 LOGICAL      :: reading_constraints=.FALSE.          ! Determine whether the reader is reading constraints at that moment
 LOGICAL      :: include_dpnd                         ! Whether to include a constraint in dependent list or not
 CHARACTER*20 :: label(1:3)

 CHARACTER*40 :: input(5) , keyword1
 CHARACTER*40 :: input_comp
 INTEGER      :: posstart,posend  ! To search the key word label
 INTEGER      :: poscomment            ! position of comment flag 
 INTEGER      :: pos_comp_strt,pos_comp_end    ! needed to analyse the component info
 INTEGER      :: ios=0, line=0, ilabel=0, i_filled_line=0   ! i_filled_line = index of non-empty line which is read by the program
 INTEGER      :: posend_keyword1,n_inp_fmt_label 
 INTEGER      :: icomp, isite_type
 INTEGER      :: icons                                      ! To count the number of constraint read 
 INTEGER      :: ninput
 LOGICAL      :: component_mismatch_flag =.FALSE.           ! TRUE when a mismatch of component is found
 LOGICAL      :: site_type_mismatch_flag =.FALSE.           ! TRUE when a mismatch of site type is found
 LOGICAL      :: ws_comp_list(1:ncomp)                      ! component list of working site   
 LOGICAL      :: constraint_found(1:ncomp,1:nsite_type)     ! To check whether constraint is read form the CONSTRAINTS file
 LOGICAL      :: file_exists                                ! To check whether CONSTRAINTS file exists in the run directory
 TYPE(CONSTRAINT),ALLOCATABLE :: chk(:)
 TYPE(DEPENDENT)  tmp_dpnd                                  !checking the repeatation of dependent
 
!Declaring counter variable
 INTEGER ii, jj,kk, i_comp_read,ll,mm,nn



!     Array allocation
 if (.not.allocated(cons)) then
    allocate (cons(1:ncomp,1:nsite_type,1:max_cons))
    do ii=1,ncomp
       do jj=1,nsite_type
          do kk=1,max_cons
             ALLOCATE(cons(ii,jj,kk)%comp_list(1:ncomp))
          enddo
       enddo
    enddo
 endif
      
 if(.not.allocated(nconstraints)) then
    allocate(nconstraints(1:ncomp,1:nsite_type))
 endif

 if(.not.allocated(h_constrain)) then
    allocate(h_constrain(1:ncomp,1:nsite_type))
 endif

 if(.not.allocated(chk)) then
    ALLOCATE(chk(1:max_cons))
    do ii=1,max_cons
       allocate(chk(ii)%comp_list(1:ncomp))
    enddo
 endif

 if (.not.allocated(dependent_list)) then
    allocate(dependent_list(1:nsite_type,1:max_dependent))
 endif

 if (.not.allocated(no_of_dependent)) then
    allocate(no_of_dependent(1:nsite_type))
 endif

      
! Array initialization
 nconstraints(:,:)=0
 inp_fmt_label(:)=default_inp_fmt_label(1:5)
 label(:)=default_label(1:3)
 ws_comp_list(:)=.FALSE.
 constraint_found(:,:) = .FALSE.
 h_constrain(:,:)=h_cons_min     
 do ii=1,ncomp
    do jj=1,nsite_type
       do kk=1,max_cons
          cons(ii,jj,kk)%comp_list(:)=.TRUE.
       enddo
    enddo
 enddo

 
 no_of_dependent(:)=0

 INQUIRE(file='CONSTRAINTS',exist=file_exists)  
 if(file_exists) then
   OPEN(104,file='CONSTRAINTS',status='old')
 else
   stop_flag=.TRUE.
   write(*,'(A)') "!!ERROR!!! CONSTRAINTS file not found! Program will abort!"
   write(*,'(A)') "Please provide CONSTRAINTS file to run simulations."
   RETURN
 endif 
   
 do while (ios.eq.0)

!     initialize mismatch flags
    component_mismatch_flag=.FALSE.
    site_type_mismatch_flag=.FALSE.

    read(104,'(A)',iostat=ios) full_input_line
   
!****  write(*,*) full_input_line, ios 

    if(ios.eq.0) then

!*     write(*,*) full_input_line,ios            
            
       line=line+1
       buffer=full_input_line
       poscomment=scan(buffer,'!#')
       if (poscomment.ne.0) then
          buffer=buffer(1:poscomment-1)
       endif
       posstart=verify(buffer,'   ')
       
       if (posstart.eq.0) then
!*     print *, "WARNING!!! Skipping blank line in ADS_FREE_EN in"
!*     x  ,line
          cycle
       endif
       
       i_filled_line=i_filled_line+1

       posend_keyword1=scan(buffer(posstart:),'   ')+posstart-2
       keyword1=buffer(posstart:posend_keyword1)

       if((i_filled_line.eq.1).and.(keyword1.eq.'Input_Format').and.(.not.reading_constraints))then
          buffer=buffer(posend_keyword1+1:)
          posstart=verify(buffer,'   ')
          
          if (posstart.eq.0) then
             Print *,'No Input_Format provided!!!Taking the default one'
             go to 1003
             
          else 
             do while ((posstart.ne.0).and.(ilabel.lt.5))
                ilabel=ilabel+1
                posend=scan(buffer(posstart:),'   ')+posstart-2
                inp_fmt_label(ilabel)=buffer(posstart:posend)
                buffer=buffer(posend+1:)
                posstart=verify(buffer,'   ') 
                if(inp_fmt_label(1).eq.'Default') then
                   inp_fmt_label(1:5)=default_inp_fmt_label(1:5)
                   go to 1003
                endif
             enddo
             n_inp_fmt_label=ilabel
             
             do ii=1,5
                if(inp_fmt_label(ii).ne.default_inp_fmt_label(ii)) then
                   if (n_inp_fmt_label.ne.5) then
                      Print *, 'You have to specify all 5 labels explicitely if you do not keep the default input format'
                      stop_flag=.TRUE.
                      exit
                   endif
                endif
             enddo
          endif
               
 1003     continue 
               
               
       else if(keyword1.eq.'Start_Constraints') then
          if(reading_constraints) then
             stop_flag=.TRUE.
             write(*,'(A)') 'ERROR in reading CONSTRAINTS at line '//trim(int2str(line))//'. Expects an End_Constraints flag before'
          else
             buffer=buffer(posend_keyword1+1:)
             posstart=verify(buffer,'   ')
             if(posstart.eq.0) then
                stop_flag=.TRUE.
                print *, 'Need the working component,site and constrained adsorption energy (optional)'
             else
                ilabel=0
                do while ((posstart.ne.0).and.(ilabel.lt.3)) 
                   ilabel=ilabel+1
                   posend=scan(buffer(posstart:),'   ')+posstart-2
                   input(ilabel)=buffer(posstart:posend)
                   buffer=buffer(posend+1:)
                   posstart=verify(buffer,'   ')
                enddo
             endif 

             ninput=ilabel 

             if(ninput.lt.2) then
                print *, 'Minimum 2 input needed in Start_Constraints (Component/Site)'
                stop_flag=.TRUE.
             endif
             
             do ilabel=1,ninput
                
                select case (label(ilabel))
                   
                case ('Component')
                   
                   if (integer_test(input(ilabel))) then
                      read(input(ilabel),*,iostat=ios) ii ! if input is an integer
                      if((ii.gt.0).and.(ii.le.ncomp)) then
                         ws_comp_list(ii)=.TRUE.
                         component_mismatch_flag=.FALSE.
!****  write(*,*) ii
                      else
                         component_mismatch_flag=.TRUE.
                      endif
                      
                   else    ! input is a string
                      
                      if(input(ilabel).eq.'All') then            
                         ws_comp_list(1:ncomp)=.TRUE.
                         
                      else !input is not all
                         comp_buff=input(ilabel)
                         pos_comp_strt=1
                         i_comp_read=0
                         do while (pos_comp_strt.ne.0) 
                            i_comp_read=i_comp_read+1        
                            pos_comp_end= scan(comp_buff(pos_comp_strt:),'| ')+pos_comp_strt-2
                            
                            input_comp=comp_buff(pos_comp_strt:pos_comp_end) 
                            
                            if(integer_test(input_comp)) then
                               read(input_comp,*,iostat=ios) ii ! if input is an integer
                               if((ii.gt.0).and.(ii.le.ncomp)) then
                                  ws_comp_list(ii)=.TRUE.
                                  component_mismatch_flag= component_mismatch_flag.or..FALSE.
!****  write(*,*) ii
                               else
                                  component_mismatch_flag=component_mismatch_flag.or..TRUE.
                               endif
                               
                            else !if input is not an integer   
                               
                               do ii=1,ncomp
                                  if(component(ii).eq.input_comp) then
                                     ws_comp_list(ii)=.TRUE.
                                     component_mismatch_flag=component_mismatch_flag.or..FALSE.
!****  write(*,*) ii
                                     exit
                                  endif
                               enddo

                               if(ii.gt.ncomp) then   
                                 component_mismatch_flag=component_mismatch_flag.or..TRUE.
                               endif                               
                            endif
                            
                            comp_buff=comp_buff(pos_comp_end+1:)
                            pos_comp_strt=verify(comp_buff,'| ') 
                            
                         enddo    
                      endif    
                   endif   
                   
                Case('Site') 
                   if(integer_test(input(ilabel))) then
                      read(input(ilabel),*,iostat=ios) ii
                      if((ii.gt.0).and.(ii.le.nsite_type)) then
                         isite_type=ii
                         site_type_mismatch_flag=.FALSE.
                      else
                         site_type_mismatch_flag=.TRUE.
                      endif
                      
                   else
                      do ii=1,nsite_type
                         if(site_type(ii).eq.input(ilabel)) then
                            isite_type=ii
                            site_type_mismatch_flag=.FALSE.
!****  write(*,*) isite_type
                            exit
                         else
                            site_type_mismatch_flag=.TRUE.
                         endif
                      enddo
                   endif
                   
                Case('Cons_En')
                   read(input(ilabel),*,iostat=ios) cons_en_temp 
                   
                Case default
                   write(*,*) 'Label not recognized in Start_Constraints block in file CONSTRAINTS in line '//trim(int2str(line))
                End Select
                
             Enddo
             
             reading_constraints=.TRUE.
             icons = 0

          endif
          
          

          
       else if (reading_constraints.and.keyword1.ne.'End_Constraints') then 
          
               
!**** 	  do icomp=1,ncomp
!****	  if (ws_comp_list(icomp)) then
!****	  write(*,*) 'icomp=',icomp,'isite_type=',isite_type
!****	  write(*,*) ws_comp_list(icomp)
!****	  endif
!****	  enddo
!****  write(*,*) 'reading_constraint loop reached'
               
          icons=icons+1 
          ilabel=0

          if(icons.gt.max_cons) then                           !checking max_cons limit
            write(*,'(A)') , 'ERROR in CONSTRAINTS!!! in line '//trim(int2str(line))
            Print *, 'More number of constraints than Max_Cons value! Increase Max_Cons!'
            stop_flag=.TRUE.
          endif

          do while ((posstart.ne.0).and.(ilabel.lt.5))
             ilabel=ilabel+1
             posend=scan(buffer(posstart:),'   ')+posstart-2
             input(ilabel)=buffer(posstart:posend)
             buffer=buffer(posend+1:)
             posstart=verify(buffer,'   ')
          enddo
          ninput=ilabel
          
          if(ninput.lt.3) then
             Print *,'Minimum 3 input (Site,X,Y) needed'
             stop_flag=.TRUE.
          endif
          
          Do ii=1,ncomp
!****	write(*,*) ws_comp_list(ii)
             if(ws_comp_list(ii)) then
              h_constrain(ii,isite_type)=cons_en_temp
                   
!****	write(*,*) ii,isite_type,h_constrain(ii,isite_type)
                   
              do ilabel=1,ninput
                 Select case (inp_fmt_label(ilabel))
                    
                 Case('Site') 
                    if(integer_test(input(ilabel))) then
!****  write(*,*) integer_test(input(ilabel)) 
                       read(input(ilabel),*,iostat=ios) jj
                       if((ii.gt.0).and.(ii.le.nsite_type)) then
                          cons(ii,isite_type,icons)%site=jj 
                          site_type_mismatch_flag=.FALSE.
                          
                          chk(icons)%site=jj
!****  write(*,*) isite_type
                       else
                          site_type_mismatch_flag=.TRUE.
                       endif
                       
                    else
                       do jj=1,nsite_type
                          if(site_type(jj).eq.input(ilabel)) then
                             cons(ii,isite_type,icons)%site=jj
                             site_type_mismatch_flag=.FALSE.

                             chk(icons)%site=jj
!****  write(*,*) isite_type
                             exit
                          else
                             site_type_mismatch_flag=.TRUE.
                          endif
                       enddo
                    endif
                    
                 Case('X')  
                    read(input(ilabel),*,iostat=ios)  cons(ii,isite_type,icons)%x 
                    
                    read(input(ilabel),*,iostat=ios)  chk(icons)%x
                    
                 Case('Y')
                    read(input(ilabel),*,iostat=ios)  cons(ii,isite_type,icons)%y
                    
                    read(input(ilabel),*,iostat=ios)  chk(icons)%y
                 Case('Component')
                    
! If component label found then initialize this particulr cons%comp_list with FALSE
                    do jj=1,ncomp 
                       cons(ii,isite_type,icons)%comp_list(jj)=.FALSE.
                    enddo
                    
                    if (integer_test(input(ilabel))) then
                       read(input(ilabel),*,iostat=ios) jj ! if input is an integer
                       if((jj.gt.0).and.(jj.le.ncomp)) then
                          cons(ii,isite_type,icons)%comp_list(jj)=.TRUE.
                          chk(icons)%comp_list(jj)=.TRUE.
                          component_mismatch_flag=.FALSE.
!****  write(*,*) ii
                       else
                          component_mismatch_flag=.TRUE.
                       endif
                       
                    else   ! input is a string
                       
                       if(input(ilabel).eq.'All') then
                          cons(ii,isite_type,icons)%comp_list(1:ncomp)=.TRUE.
                          chk(icons)%comp_list(1:ncomp)=.TRUE.
                          
                       else !input is not all
                          comp_buff=input(ilabel)
                          pos_comp_strt=1
                          i_comp_read=0
                          do while (pos_comp_strt.ne.0)
                             i_comp_read=i_comp_read+1
                             pos_comp_end=scan(comp_buff(pos_comp_strt:),'| ')+pos_comp_strt-2
                             
                             input_comp=comp_buff(pos_comp_strt:pos_comp_end)
                             
                             if(integer_test(input_comp)) then
                                read(input_comp,*,iostat=ios)jj ! if input is an integer
                                if((jj.gt.0).and.(jj.le.ncomp)) then
                                   cons(ii,isite_type,icons)%comp_list(jj)=.TRUE.
                                   chk(icons)%comp_list(jj)=.TRUE.
                                   component_mismatch_flag=component_mismatch_flag.or..FALSE.
!****  write(*,*) ii
                                else
                                   component_mismatch_flag= component_mismatch_flag.or..TRUE.
                                endif
                                
                             else !if input is not an integer   
                                
                                do jj=1,ncomp
                                   if(component(jj).eq.input_comp) then
                                      cons(ii,isite_type,icons)%comp_list(jj)=.TRUE.
                                    chk(icons)%comp_list(jj)=.TRUE.
                                      component_mismatch_flag=component_mismatch_flag.or..FALSE.
!****  write(*,*) ii
                                      exit
                                   endif
                                enddo
                              
                                if (jj.gt.ncomp) then
                                  component_mismatch_flag=component_mismatch_flag.or..TRUE.
                                endif
                             endif
                             
                             comp_buff=comp_buff(pos_comp_end+1:)
                             pos_comp_strt=verify(comp_buff,'| ')
                             
                          enddo
                       endif
                    endif
                    
                 Case('Present')
                    read(input(ilabel),*,iostat=ios) cons(ii,isite_type,icons)%switch   
                    read(input(ilabel),*,iostat=ios)chk(icons)%switch
                 Case default
                    Print *, 'Keyword in Input_Format not recognised'
                    stop_flag=.TRUE.
                 End select
                 
              Enddo
           Endif
        Enddo
! Checking repeatation error in the input file             
        do ii=1,icons-1
           if((chk(icons)%site.eq.chk(ii)%site).and.(chk(icons)%x.eq.chk(ii)%x).and.(chk(icons)%y.eq.chk(ii)%y)) then
                if(chk(icons)%switch.eqv.chk(ii)%switch) then
              
                  write(*,'(A)') 'ERROR in CONSTRAINTS!!! at line '//trim(int2str(line))
                  write(*,'(A)') 'Repeatation of constraint not allowed! All of "Site" "X" and "Y" and "Switch" & 
               &  can not repeat in a same Constraint block'
                  stop_flag=.TRUE.
           
               else 
                   do jj=1,ncomp
                      if(chk(icons)%comp_list(jj).eqv.chk(ii)%comp_list(jj)) then
                        write(*,'(A)') 'ERROR in CONSTRAINTS!!! at line '//trim(int2str(line))
                        write(*,'(A)') 'Repeatation of constraint not allowed! As all "Site" "X" "Y" are same, & 
                    &   Same component cannot appear in both line'
                        stop_flag=.TRUE.
                        exit
                      endif
                    enddo
               endif 
           endif
        enddo
 
!     iconstraint_read=iconstraint_read+1
        reading_constraints=.TRUE.
        
     else if (keyword1.eq.'End_Constraints') then
        if(reading_constraints) then
           do ii=1,ncomp
              if(ws_comp_list(ii)) then
                 nconstraints(ii,isite_type)=icons
                 constraint_found(ii,isite_type)=.TRUE.
              endif
           enddo
           reading_constraints=.FALSE.
           ws_comp_list(:)=.FALSE.
           cons_en_temp=0.0d0
        else
           stop_flag=.TRUE.
           write(*,'(A)') 'ERROR in CONSTRAINTS!!! at line '//trim(int2str(line))//'. Expects an Start_Constraints flag before'
        endif
       
   
     else
        stop_flag=.TRUE.
        write(*,'(A)') 'ERROR in CONSTRAINTS!!! at line'//trim(int2str(line))//'. Invalid keyword'
        
        
     endif
     
     if(component_mismatch_flag) then
        write(*,'(A)')'ERROR in CONSTRAINTS!!! Component mismatch found in CONSTRAINTS at line '//trim(int2str(line))
        stop_flag=.TRUE.
     endif
     
     if(site_type_mismatch_flag) then
        write(*,'(A)') 'ERROR in CONSTRAINTS!!! Site type mismatch found in CONSTRAINTS at line '//trim(int2str(line))
        stop_flag=.TRUE.
     endif
     
  endif
 enddo 
      
! Checking consistency between ADS_FREE_EN and CONSTRAINTS file after completion of reading
     
 deallocate(chk)

 do ii=1,ncomp
    do jj=1,nsite_type 

       if(h_cons_min.gt.h_constrain(ii,jj)) then
         h_cons_min=h_constrain(ii,jj)
       endif

       if(constraint_present(ii,jj).neqv.constraint_found(ii,jj)) then

          write(*,'(A)') 'WARNING!!!Constraint mismatch found between ADS_EN and CONSTRAINTS'
          write(*,'(A)') 'Component= '//trim(component(ii))//' Site= '//trim(site_type(jj))

          Print *, 'Constraint present in ADS_FREE_EN?= ',constraint_present(ii,jj)
          Print *, 'Constraint found in CONSTRAINTS?= ',constraint_found(ii,jj)
          if(constraint_present(ii,jj)) then
            Print *, 'ERROR in CONSTRAINTS!!!Constraint not found! Will abort'
            stop_flag=.TRUE.
          else
            Print *, 'Continuing neglecting the constraint data provided in CONSTRAINTS'
          endif
       endif
    enddo
 enddo

     
! Making dependent's list

 if(.not.stop_flag) then
 do ii=1,ncomp
    do jj=1,nsite_type
       do kk=1,nconstraints(ii,jj)
          ll=cons(ii,jj,kk)%site
          tmp_dpnd%site = jj
          tmp_dpnd%x    = -1*cons(ii,jj,kk)%x
          tmp_dpnd%y    = -1*cons(ii,jj,kk)%y
          include_dpnd=.TRUE.

          if(no_of_dependent(ll).ge.1) then
           do mm=1,no_of_dependent(ll)
             if((dependent_list(ll,mm)%x.eq.tmp_dpnd%x).and.(dependent_list(ll,mm)%y.eq.tmp_dpnd%y).and. &
            &  (dependent_list(ll,mm)%site.eq.tmp_dpnd%site)) then 
                include_dpnd=.FALSE.
                exit
             endif 
           enddo
          endif      

         if(include_dpnd) then
         no_of_dependent(ll)=no_of_dependent(ll)+1
          if(no_of_dependent(ll).gt.max_dependent) then
           write(*,*) 'ERROR in Input!!! Please increase the Max_Dependent value in SIMULATION_CONTROL'
           stop_flag=.TRUE.
           exit
          endif
         dependent_list(ll,no_of_dependent(ll))%x=tmp_dpnd%x
         dependent_list(ll,no_of_dependent(ll))%y=tmp_dpnd%y
         dependent_list(ll,no_of_dependent(ll))%site=tmp_dpnd%site
         endif

       enddo
    enddo
 enddo
 endif

!check debug
!*      do ii=1,nsite_type
!*       write(*,*) 'number of dependent sites of'//site_type(ii),'are',
!*     x no_of_dependent(ii)
!*       if(no_of_dependent(ii).gt.0) then
!*        write(*,*) 'X    Y     Site'
!*        do jj=1,no_of_dependent(ii)
!*           write(*,*) dependent_list(ii,jj)%x,dependent_list(ii,jj)%y,
!*     x  site_type(dependent_list(ii,jj)%site)
!*        enddo
!*       endif
!*      enddo
!************************************************************************************
        
          
             

               
 RETURN 
 END SUBROUTINE read_constraints

      
      
      
                
