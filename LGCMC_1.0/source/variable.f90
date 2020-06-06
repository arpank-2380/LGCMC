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


 MODULE variable    !All global variables
 USE derived_data
 IMPLICIT NONE
 SAVE
!     Simulation control related variable     
 DOUBLE PRECISION ::  t !Temperature
 INTEGER          :: neqbmstep=1000, ntotalstep
 INTEGER(KIND=8)  :: nsimulstep=1000
 INTEGER          :: nxsite=6, nysite=10     !initialized for Mg-MOF system ! simulation of single pore
 INTEGER(KIND=8)  :: nsite   !Total number of a particular type of site
 INTEGER          :: ncomp=1   !number of component default pure system
 INTEGER          :: nsite_type =1  !number of different sites types
!      DOUBLE PRECISION, PARAMETER :: x1=0.1 !Mole fraction of component 1
!      DOUBLE PRECISION, PARAMETER :: x(2)=(/x1,1-x1/) !Mole fractions of components in gas
 LOGICAL          :: restart=.FALSE. !Do you want to restart from some saved configuration       
 LOGICAL          :: corr_func=.TRUE. !Do you like to calculate the default correlation functions

 INTEGER          :: corx1,corx2,cory1,cory2 ! pair correlation lattice distance cut off

 CHARACTER        :: eos*30='Ideal'     !Equation of state being used to calculate fugacity (chemical potential)
 CHARACTER*20,ALLOCATABLE  :: component(:),site_type(:)  !Reading the components name and site types in an array
 TYPE (COMP_INFO), ALLOCATABLE :: comp_information(:)    !Critical temperature, pressure and accentric factor
 LOGICAL          ::  gas_int=.FALSE.                    !If .TRUE. user can give input for int_param_a and int_param_b through nput file GAS_INT 
 DOUBLE PRECISION, ALLOCATABLE :: int_param_a(:,:)       !Gas phase binary interaction parameters for mixed term "A"
 DOUBLE PRECISION, ALLOCATABLE :: int_param_b(:,:)       !Gas phase binary interaction parameters for mixed term "B"


 LOGICAL          :: stop_flag=.FALSE.                   !program stop flag
 LOGICAL          :: displacement=.TRUE.                 ! Whether displacement routine is on in simulation
 INTEGER          :: disp_start                          ! When diplsce ment move will be started
 INTEGER          :: x_disp=3                            ! displacement move upto +- x_disp lattice distance
 INTEGER          :: y_disp=3                            ! displacement move upto +- y_disp lattice distance
 INTEGER          :: disp_box                            ! Number of exchange possible in the displacement box
 INTEGER          :: maxtrial=3                          ! Maximum # of trial move= maxtrial*4*x_disp*y_disp*nsite_type
 INTEGER          :: tdump=10                            !dump step to calculate average
 INTEGER          :: maxdisp                             ! Maximum number of displace attempt
 INTEGER          :: iseed=8745                          !integer seed to generate random number

 DOUBLE PRECISION, ALLOCATABLE :: h(:,:)                 !adsorption free energy

!     Constraint specific variables
 LOGICAL,ALLOCATABLE :: constraint_present(:,:)            ! site and component specific constraints in adsorption energy
 LOGICAL             :: constraint_file_needed=.FALSE.     !  Determine whether constraint_file is needed or not             
 INTEGER, ALLOCATABLE :: nconstraints(:,:)                 ! Total number of constraints
 INTEGER :: max_cons=10                                       ! maximum number of constraints
 TYPE (CONSTRAINT), ALLOCATABLE :: cons(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: h_constrain(:,:)
 DOUBLE PRECISION :: h_cons_min
 TYPE(DEPENDENT), ALLOCATABLE :: dependent_list(:,:)       !List of dependent sites
 INTEGER, ALLOCATABLE :: no_of_dependent(:)                ! no of dependent site of a particular site.
 INTEGER :: max_dependent=10                               ! Maximum number of allowed dependent
 TYPE(DEPENDENT), ALLOCATABLE :: virtual_sites(:)          !list of virtual sites
 
!     Lateral interaction specific variable
 LOGICAL :: lat_int=.TRUE.
 INTEGER :: xcut=1,ycut=1                             ! cut off distance along two axis.
 DOUBLE PRECISION, ALLOCATABLE :: j(:,:,:,:,:,:)        ! lateral interaction parameter
 DOUBLE PRECISION, ALLOCATABLE :: jb(:,:,:,:,:,:)       ! lateral interaction for lattice assymmetry
 LOGICAL,ALLOCATABLE :: read_lat_int(:,:,:,:,:,:) ! Whether the lateral interaction array element has been read
 LOGICAL,ALLOCATABLE :: asym(:,:,:,:,:,:)         ! Lattice_assymmetry present or not

!     Composition of a particular site
 INTEGER, ALLOCATABLE :: confg(:,:,:) 

 END MODULE variable 
