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


 SUBROUTINE read_gas_int(gas_param_a,filename,unit_no,file_exists)
 IMPLICIT NONE
 LOGICAL, INTENT (IN) ::  gas_param_a
 CHARACTER*9, INTENT(IN) ::  filename
 INTEGER, INTENT(IN)      ::  unit_no
 LOGICAL, INTENT(OUT)     ::  file_exists
 CHARACTER*100            ::  full_input_line, buffer
 CHARACTER*40             ::  input(2)                                          ! Temporary variables to read the input components
 DOUBLE PRECISION         ::  input_value=0.0d0                                       ! Temporary variable to read the  input data
 INTEGER                  ::  ios=0, ios2=0, line=0                             ! input-output status specifier and line counter 
 INTEGER                  ::  ii,jj, icounter,posstart,poscomment,icomp(2)
 LOGICAL                  ::  component_mismatch(2)=(/.FALSE.,.FALSE./) 
 LOGICAL                  ::  component_mismatch_flag=.FALSE.                          

 
 INQUIRE(file=filename,exist=file_exists)
 if (file_exists) then
    OPEN(unit=unit_no,file=filename,status='old')
 else
    RETURN
 endif

! Initializing the variables again
 ios=0
 line=0
 ios2=0
 component_mismatch(:)=(/.FALSE.,.FALSE./)
 component_mismatch_flag=.FALSE.

 do while (ios.eq.0)
    read(unit_no,'(A)',iostat=ios) full_input_line
    if (ios.eq.0) then
       line=line+1
       buffer=full_input_line
       poscomment=scan(buffer,'!#')
       if (poscomment.ne.0) then
          buffer=buffer(1:poscomment-1)
       endif
       posstart=verify(buffer,'   ')
 
       if(posstart.eq.0) then
         cycle
       else    
         read(buffer,*,iostat=ios2) input(1), input(2), input_value
       endif

       if (ios2.eq.0) then
          Do icounter=1,2
             if (integer_test(input(icounter))) then
                read(input(icounter),*,iostat=ios2) ii
                if((ii.gt.0).and.(ii.le.ncomp)) then
                    icomp(icounter)=ii
                    component_mismatch(icounter)=.FALSE.
                else
                    component_mismatch(icounter)=.TRUE.
                endif
             else
                do ii=1,ncomp
                   if(component(ii).eq.input(icounter)) then
                      icomp(icounter)=ii
                      component_mismatch(icounter)=.FALSE.
                      exit
                   else
                      component_mismatch(icounter)=.TRUE.
                   endif
                enddo
             endif
          Enddo

          component_mismatch_flag=component_mismatch(1).or.component_mismatch(2)  
         
          if (.not.component_mismatch_flag) then
             if (gas_param_a) then
                 int_param_a(icomp(1),icomp(2)) = input_value
                 int_param_a(icomp(2),icomp(1)) = input_value
             else
                 int_param_b(icomp(1),icomp(2)) = input_value
                 int_param_b(icomp(2),icomp(1)) = input_value
             endif
          else
             write(*,*)
             write(*,'(A)') "!!! ERROR in "//trim(filename)//" !!! Component type mismatch occured at line "//trim(int2str(line)) 
             stop_flag=.TRUE.
             RETURN
          endif

       else
           write(*,'(A)') "!!! ERROR reading "//trim(filename)//" !!! See below for suggestions."
           write(*,'(A)') "Usage: Component1{<Name:String> or <Index: Integer>} Component2{<Name:String> or <Index: Integer>} Interaction_Parameter{<FLOAT>}"
           stop_flag=.TRUE. 
           RETURN
       endif
    endif
enddo

RETURN
END SUBROUTINE read_gas_int     

