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


 MODULE energy_calculation
 USE variable, ONLY : lat_int
 USE calc_ads_en, ONLY : ads_en
 USE calc_lat_int_en, ONLY : lat_int_en

 CONTAINS
 FUNCTION energy(ixsite,iysite,isite_type,icomp)
 DOUBLE PRECISION energy
 INTEGER, INTENT(IN) :: ixsite,iysite,isite_type,icomp
 if(lat_int) then
    energy=ads_en(ixsite,iysite,isite_type,icomp)+lat_int_en(ixsite,iysite,isite_type,icomp)
 else
    energy=ads_en(ixsite,iysite,isite_type,icomp)
 endif
 RETURN
 END FUNCTION energy

 END MODULE energy_calculation
      
