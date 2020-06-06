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


MODULE logo_module

CONTAINS

SUBROUTINE credit_print(label)
IMPLICIT NONE
INTEGER, INTENT(IN) :: label
write(label,'(59A)') "###########################################################"
write(label,'(59A)') "#              This output is generated by                #"
write(label,'(59A)') "#                ╦    ╔═╗  ╔═╗  ╔╦╗  ╔═╗                  #" 
write(label,'(59A)') "#                ║    ║ ╦  ║    ║║║  ║                    #"
write(label,'(59A)') "#                ╩═╝  ╚═╝  ╚═╝  ╩ ╩  ╚═╝                  #"
write(label,'(59A)') "#       A Lattice Grand Canonical Monte Carlo Code        #" 
write(label,'(59A)') "#                Written by Arpan Kundu                   #"
write(label,'(59A)') "#                        at                               #"
write(label,'(59A)') "#            Humboldt University of Berlin                #"
write(label,'(59A)') "###########################################################"
END SUBROUTINE credit_print

END MODULE logo_module
