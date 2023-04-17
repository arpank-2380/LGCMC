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


 MODULE read_input_files
 
 USE variable 
!****      LOGICAL,ALLOCATABLE :: read_lat_int(:,:,:,:,:,:) ! Whether the lateral interaction array element has been read   
 DOUBLE PRECISION, ALLOCATABLE :: j_raw(:,:,:,:,:,:)  !Dimerization energy supplied by user
!     Gas related data arrays
 DOUBLE PRECISION, ALLOCATABLE :: gas_p(:)
 DOUBLE PRECISION, ALLOCATABLE :: gas_comp(:,:)
 DOUBLE PRECISION, ALLOCATABLE :: gas_fug(:,:)
!     gas.dat related variable
 INTEGER :: npoint                                    !Number of point given in gas.dat 
 LOGICAL :: gas_int_a_exists=.FALSE., gas_int_b_exists=.FALSE.  !.TRUE. if GAS_INT file or GAS_INT_B file exists          
 LOGICAL :: eos_implemented = .FALSE.                 ! Checks whether the asked EOS is implemented in the code 

 CONTAINS
  
   INCLUDE 'inp_read_routines/read_simulation_control.f90'

   INCLUDE 'inp_read_routines/read_ads_free_en.f90'

   INCLUDE 'inp_read_routines/read_constraints.f90'

   INCLUDE 'inp_read_routines/read_lat_int.f90'

   INCLUDE 'inp_read_routines/read_gas_data.f90'

   INCLUDE 'inp_read_routines/read_gas_int.f90'

 FUNCTION integer_test(string)
 LOGICAL :: integer_test
 CHARACTER*40 string
 if(verify(string,'0123456789 ').eq.0) then
   integer_test=.TRUE.
 else
   integer_test=.FALSE.
 endif
 RETURN
 END FUNCTION integer_test

 FUNCTION int2str(input_int)
 CHARACTER*20 int2str
 INTEGER, INTENT(IN) :: input_int
 INTEGER :: ios
 write(int2str,*,iostat=ios) input_int
 if (ios.ne.0) then
    write(*,*) 'Error reading input/Input is not an integer' 
 endif
 call remove_leading_blanks(int2str)
 RETURN   
 END FUNCTION int2str 

 FUNCTION dble2str(input_dble)
 CHARACTER*30 dble2str
 DOUBLE PRECISION, INTENT(IN) :: input_dble
 INTEGER :: ios
 write(dble2str,'(G20.10)',iostat=ios) input_dble
 if(ios.ne.0) then
    write(*,*) 'Error reading input/ Input is no a DBLE number'
 endif
 call remove_leading_blanks(dble2str)
 RETURN
 END FUNCTION dble2str

 SUBROUTINE remove_leading_blanks(input_string)

 INTEGER i,j,strln
 CHARACTER (*) input_string

 if(input_string(1:1).ne.' ') then
    return
 endif

 strln=len(input_string)

 do i=1,strln
   if(input_string(i:i).ne.' ') then
     exit
   endif
 enddo

 input_string(1:strln-i+1)=input_string(i:strln)

 do j=strln-i+2,strln
    input_string(j:j)=' '
 enddo

 RETURN    
 END SUBROUTINE remove_leading_blanks

 END MODULE read_input_files 

      

