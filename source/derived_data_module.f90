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


 MODULE derived_data
!     Declaration of new data types
 TYPE CONSTRAINT
 SEQUENCE
     INTEGER :: site
     INTEGER :: x,y
     LOGICAL :: switch=.TRUE.
     LOGICAL, ALLOCATABLE :: comp_list(:)
 ENDTYPE CONSTRAINT

 TYPE DEPENDENT
 SEQUENCE
     INTEGER :: x,y,site
 ENDTYPE DEPENDENT

 TYPE COMP_INFO
 SEQUENCE
      DOUBLE PRECISION :: pc,tc,acc_fac
 ENDTYPE COMP_INFO

 END MODULE derived_data
