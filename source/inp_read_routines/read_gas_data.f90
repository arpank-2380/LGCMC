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


 SUBROUTINE read_gas_data
 IMPLICIT NONE

! Input related variables
 CHARACTER*300  :: input_string,buffer
 CHARACTER*40,ALLOCATABLE   :: label1(:), label2(:), input(:)
 CHARACTER*40   :: temp_label
 INTEGER, ALLOCATABLE :: label_comp(:)                        !To store the component part of the label 
 INTEGER        :: posstart,posend                            ! To search the key word label
 INTEGER        :: poscomment                                 ! position of comment flag 
 INTEGER        :: temp_label_len, delimeter_pos              ! to break keyword y_ and Phi_
 INTEGER        :: ios=0                                      ! input output status                 
 INTEGER        :: ilabel=0, nlabel                           ! total number of label(key-word) 
 INTEGER        :: icomp=1,jcomp=1, temp_comp                 ! index of component and a temporary variable to store component-index
 INTEGER        :: maxlabel                                   ! maximum number of allowed labels
 INTEGER        :: label_p=0, label_y=0, label_phi=0          ! Number of supplied labels to check they are enough
 INTEGER        :: y2calc=1                                   ! y of which component-index must be calculated and not supplied
!Temporary storage variables

 LOGICAL        :: component_mismatch_flag =.FALSE.           ! TRUE when a mismatch of component is found
 LOGICAL        :: file_exists                                !TRUE when ADS_FREE_EN file exists
 LOGICAL        :: data_stop_flag=.FALSE.                     !Column-data related stop flag
 LOGICAL,ALLOCATABLE :: comp2calc(:)                          !TRUE when the composition is supposed to be calculated
      
!Declaring counter variables
 INTEGER ii,ipoint,ninput     
 INTEGER read_p, read_comp, read_phi                          !counter variable to chek how many data are read from each line
!Printing advance option
 DOUBLE PRECISION :: y_sum                                    ! variable to store sum of composition of all components

 maxlabel=2*ncomp

!Opening files
  
 INQUIRE(file='GAS',exist=file_exists)
 if (file_exists) then
    Open(105,file='GAS',status='old')
 else
    stop_flag=.TRUE.
    write(*,'(A)') "!!!ERROR!!! GAS file not found! Program will abort!"
    write(*,'(A)') "Please provide GAS file to continue with the simulations."
    RETURN
 endif
 
 ALLOCATE(label1(1:maxlabel),label2(1:maxlabel),input(1:maxlabel))
 ALLOCATE(label_comp(1:maxlabel))
 ALLOCATE(comp2calc(1:ncomp))
! ALLOCATE(int_param_a(1:ncomp,1:ncomp),int_param_b(1:ncomp,1:ncomp))

 !Initializing the interaction parameter to be 0
! do icomp=1,ncomp
!    do jcomp=1,ncomp
!       int_param_a(icomp,jcomp)=0.0d0
!       int_param_b(icomp,jcomp)=0.0d0
!    enddo
! enddo 
 
  read(105,'(A)',iostat=ios) input_string

  if(ios.eq.0) then
    buffer=input_string
    poscomment=scan(buffer,'!#')
    if (poscomment.ne.0) then
       buffer=buffer(1:poscomment-1)
    endif 
    posstart=verify(buffer,'   ')   

    if (posstart.eq.0) then
      stop_flag=.TRUE.
      write(*,'(A)') "ERROR!!! 1st line in GAS is blank. This line is reserved for keywords!"         
      RETURN
    endif
  else 
    write(*,'(A)')"ERROR reading keywords in GAS"
  endif


!Reading labels (key-words)

 do while(posstart.ne.0)
    ilabel=ilabel+1
    if(ilabel.gt.maxlabel) then
      stop_flag=.TRUE.
      write(*,'(A)') "ERROR in GAS!!! Maximum number("//trim(int2str(maxlabel))//") of allowed keywords exceeded."
      RETURN
    else
      posend=scan(buffer(posstart:),'   ')+posstart-2
      temp_label=buffer(posstart:posend)
      buffer=buffer(posend+1:)
      posstart=verify(buffer,'   ') 
      temp_label_len=len_trim(temp_label)
      delimeter_pos=scan(temp_label,'_',.TRUE.)
      if(delimeter_pos.ne.0) then
        label1(ilabel)=temp_label(1:delimeter_pos-1)
        label2(ilabel)=temp_label(delimeter_pos+1:temp_label_len)
      else
        label1(ilabel)=temp_label
        label2(ilabel)=''
      endif
    endif
 enddo
    
 nlabel=ilabel 


!  Initializing comp_to_calc
 do icomp=1,ncomp
    comp2calc(icomp)=.TRUE.
 enddo
      
!  Checking supplied labels
 do ilabel=1,nlabel
    select case(label1(ilabel))
    case('P') 
                label_p=label_p+1
                label_comp(ilabel)=0
    case('y') 
                label_y=label_y+1
                if(integer_test(label2(ilabel))) then     !if componet index is given

                  read(label2(ilabel),*) temp_comp
                  if((temp_comp.gt.0).and.(temp_comp.le.ncomp))then
                     label_comp(ilabel)=temp_comp
                     comp2calc(temp_comp)=.FALSE.
                  else
                     component_mismatch_flag=.TRUE.
                     write(*,*) "ERROR in GAS!!! Component type mismatch found at "//trim(int2str(ilabel))//"-th keyword."
                  endif

                else
                  do ii=1,ncomp
                     if(component(ii).eq.label2(ilabel)) then
                        label_comp(ilabel)=ii
                        comp2calc(ii)=.FALSE.
                        exit
                     endif
                  enddo
                  if(ii.gt.ncomp) then   
                    component_mismatch_flag=.TRUE.
                    write(*,*) "ERROR in GAS!!! Component type mismatch found at "//trim(int2str(ilabel))//"-th keyword."
                  endif
                endif

    case('Phi') 
                label_phi=label_phi+1
                if(integer_test(label2(ilabel))) then
                  read(label2(ilabel),*) temp_comp
                  if((temp_comp.gt.0).and.(temp_comp.le.ncomp))then
                     label_comp(ilabel)=temp_comp
                  else
                     component_mismatch_flag=.TRUE.
                     write(*,*) "ERROR in GAS!!! Component type mismatch found at "//trim(int2str(ilabel))//"-th keyword."
                  endif

                else
                  do ii=1,ncomp
                     if(component(ii).eq.label2(ilabel)) then
                        label_comp(ilabel)=ii
                        exit
                     endif
                  enddo
                  if(ii.gt.ncomp) then
                    component_mismatch_flag=.TRUE.
                    write(*,*) "ERROR in GAS!!! Component type mismatch found at "//trim(int2str(ilabel))//"-th keyword."
                  endif
                endif

    case default
         stop_flag=.TRUE.
         write(*,'(A)')"ERROR in GAS!!! At 1st line "//trim(int2str(ilabel))//"-th keyword is INVALID!"
         RETURN    
    end select
 enddo

 if(component_mismatch_flag) then
   stop_flag=.TRUE.
   RETURN
 endif

 if(eos_implemented) then
   write(*,'(A)') "Expected data columns in GAS: total pressure (Keyword: P) and " & 
 &  //trim(int2str(ncomp-1))//" composition(s) (Keyword: y_)."
   if(label_phi.gt.0) then
     write(*,'(A)') "WARNING! You have chosen an implemented EOS according to which fugacity coefficients will be calculated. Supplied fugacity coefficients will not be used."
   endif
 else
   write(*,'(A)') "Expected data columns in GAS: total pressure (Keyword: P) and "//trim(int2str(ncomp-1))//&
   &" composition(s) (Keyword: y_) and "//trim(int2str(ncomp))//" fugacity coefficient(s) (Keyword: Phi_)."
   if(label_phi.ne.ncomp) then
     stop_flag=.TRUE.
     write(*,'(A)') "ERROR in GAS!!! Keyword Phi_appeared "//trim(int2str(label_phi))//" number of times." 
     write(*,'(A)') "Keyword Phi_ must appear "//trim(int2str(ncomp))//" number of times."
     write(*,'(A)') "You must provide only 1 column of fugacity coefficient for each gas."
     RETURN
   endif
 endif


 if(label_p.ne.1) then
   stop_flag=.TRUE.
   write(*,'(A)') "ERROR in GAS!!! Keyword P appeared "//trim(int2str(label_p))//" number of times." 
   write(*,'(A)') "Keyword P must appear but for only once."
   write(*,'(A)') "Only 1 column of total pressure (P) data is allowed"
   RETURN
 endif

 if(label_y.ne.(ncomp-1)) then
   stop_flag=.TRUE.
   write(*,'(A)') "ERROR in GAS!!! Keyword y_ appeared "//trim(int2str(label_y))//" number of times." 
   write(*,'(A)') "Keyword y_ must appear "//trim(int2str(ncomp-1))//" number of times."
   write(*,'(A)') "Only "//trim(int2str(ncomp-1))//" column(s) of gas composition data is(are) allowed."
   RETURN
 endif 


!     Whcih composition (y) is not supplied and must be calculated? 
 do icomp=1,ncomp
    if(comp2calc(icomp)) then
      y2calc=icomp
    endif         
 enddo 

!     Reading number of point to simulate
 read(105,*) npoint 
!     Allocation of the arrays to store gas_data
 ALLOCATE(gas_p(1:npoint),gas_comp(1:npoint,1:ncomp))
 ALLOCATE(gas_fug(1:npoint,1:ncomp))
 
!     Reading column of data
 do ipoint=1,npoint                                   !start of reading data rows

    read(105,'(A)',iostat=ios) input_string
    read_p=0                                           !initialization of the counters at the begenning of each line
    read_comp=0
    read_phi=0

    if(ios.eq.0) then
      buffer=input_string
      poscomment=scan(buffer,'!#')
      if (poscomment.ne.0) then
         buffer=buffer(1:poscomment-1)
      endif
      posstart=verify(buffer,'   ')

      if (posstart.eq.0) then
        data_stop_flag=.TRUE.
        write(*,'(A)') "ERROR in GAS!!! Blank line encountered at line "//trim(int2str(ipoint+2))//"!!!"
      endif
    else
      data_stop_flag=.TRUE.
      write(*,'(A)')"ERROR reading line "//trim(int2str(ipoint+2))//" of GAS"
    endif

    ilabel=0
    do while((ilabel.lt.nlabel).and.(posstart.ne.0))             
       ilabel=ilabel+1                                    !reading the input string
       posend=scan(buffer(posstart:),'   ')+posstart-2
       input(ilabel)=buffer(posstart:posend)
       buffer=buffer(posend+1:)
       posstart=verify(buffer,'   ')
    enddo
    ninput=ilabel
    if (ninput.lt.nlabel) then
       write(*,'(A)')'WARNING! less number of input in line '//trim(int2str(ipoint+2))//' than label specified in GAS'
    endif

    do ilabel=1,ninput
       select case (label1(ilabel))
       case('P') 
                 read(input(ilabel),*,iostat=ios) gas_p(ipoint)
                 if(ios.ne.0) then
                   data_stop_flag=.TRUE.
                   write(*,'(A)') "ERROR reading P at line "//trim(int2str(ipoint+2))//"!!! Expects a FLOAT."
                   exit
                 else
                   read_p=read_p+1
                 endif

       case('y') 
                 read(input(ilabel),*,iostat=ios) gas_comp(ipoint,label_comp(ilabel))
                 if(ios.ne.0) then
                   data_stop_flag=.TRUE.
                   write(*,'(A)') "ERROR reading y_"//trim(label2(ilabel))//"at line "//trim(int2str(ipoint+2))//"!!! Expects a FLOAT."
                   exit
                 else
                   read_comp=read_comp+1
                 endif

                 if(gas_comp(ipoint,label_comp(ilabel)).lt.0.0d0) then
                   data_stop_flag=.TRUE.
                   write(8,'(A)') "ERROR in supplied y_"//trim(label2(ilabel))//"at line "//trim(int2str(ipoint+2))// & 
                 & "!!! This value can't be negative!!!"
                   exit
                 endif
     
       case('Phi')
                  read(input(ilabel),*,iostat=ios) gas_fug(ipoint,label_comp(ilabel))
                  if(ios.ne.0) then
                    data_stop_flag=.TRUE.
                    write(*,'(A)') "ERROR reading Phi_"//trim(label2(ilabel))//"at line "//trim(int2str(ipoint+2))// & 
                   & "!!! Expects a FLOAT."
                    exit
                  else
                    read_phi=read_phi+1
                  endif
       end select
    enddo

    if (ios.eq.0) then
      y_sum=0.0d0                               !Now calculating the composition of the remaining component
      do icomp=1,ncomp
        if(icomp.ne.y2calc) then
          y_sum=y_sum+gas_comp(ipoint,icomp)
        endif
      enddo
      if(y_sum.gt.1.0d0) then
        data_stop_flag=.TRUE.
        write(*,'(A)') "ERROR in GAS!!! Sum of the supplied gas compositions is exceeding 1.0 at line "// & 
       & trim(int2str(ipoint+2))//"!!!"
      else
        gas_comp(ipoint,y2calc)=1.0d0-y_sum
      endif
    endif

! Checking if sufficient data is available in each line
    if(read_p.ne.1) then
      data_stop_flag=.TRUE.
      write(*,'(A)') "ERROR at line "//trim(int2str(ipoint+2))//" of GAS. Required number of pressure data not found."
    endif

    if(read_comp.ne.(ncomp-1)) then
      data_stop_flag=.TRUE.
      write(*,'(A)') "ERROR at line "//trim(int2str(ipoint+2))//" of GAS. Required number of composition data not found."
    endif

    if((.not.eos_implemented).and.(read_phi.ne.ncomp)) then
      data_stop_flag=.TRUE.
      write(*,'(A)') "ERROR at line "//trim(int2str(ipoint+2))//" of GAS. Required number of fugacity coefficients not found."
    endif

 end do                                       !end of reading data rows

 if(data_stop_flag) then
   stop_flag=.TRUE.
   write(*,'(A)') "ERROR encountered in column data of GAS!!!"
 endif

 DEALLOCATE(label1,label2,input,label_comp,comp2calc)

 RETURN
 END SUBROUTINE read_gas_data 
