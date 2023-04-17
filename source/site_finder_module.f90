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


! Includes the codes containing to find the site

 MODULE site_finder_module

 USE variable, ONLY : nxsite,nysite
 IMPLICIT NONE
 PUBLIC  :: site_finder     
 PRIVATE :: ixprev,ixsucc,iyprev,iysucc


!Initialize jxsite at ixsite and jysite at iysite
      
 CONTAINS

 SUBROUTINE site_finder(ixsite,iysite,ix,iy,jxsite,jysite)
    INTEGER, INTENT(IN) :: ixsite,iysite,ix,iy
    INTEGER             :: jxsite,jysite
    INTEGER             :: ii,jj                 !counter variables

    jxsite=ixsite
    jysite=iysite
    If (ix.ge.0) then
       Do ii=1,ix
          jxsite=ixsucc(jxsite)
       Enddo
    Else
       Do ii=1,Abs(ix)
          jxsite=ixprev(jxsite)
       Enddo
    Endif

    If (iy.ge.0) then
       Do ii=1,iy
          jysite=iysucc(jysite)
       Enddo
    Else
       Do ii=1,Abs(iy)
          jysite=iyprev(jysite)
       Enddo
    Endif
 END SUBROUTINE site_finder

 FUNCTION ixprev(ixsite)
 INTEGER ixprev
 INTEGER, INTENT(IN) :: ixsite
 if ((ixsite-1).ne.0) then
    ixprev=ixsite-1
 else
    ixprev=nxsite
 endif
 RETURN
 END FUNCTION ixprev

 FUNCTION ixsucc(ixsite)
 INTEGER ixsucc
 INTEGER, INTENT(IN) :: ixsite
 if((ixsite+1).le.nxsite) then
   ixsucc=ixsite+1
 else
   ixsucc=1
 endif
 RETURN
 END FUNCTION ixsucc

 FUNCTION iyprev(iysite)
 INTEGER iyprev
 INTEGER, INTENT(IN) :: iysite
 if ((iysite-1).ne.0) then
    iyprev=iysite-1
 else
    iyprev=nysite
 endif
 RETURN
 END FUNCTION iyprev

 FUNCTION iysucc(iysite)
 INTEGER iysucc
 INTEGER, INTENT(IN) :: iysite
 if ((iysite+1).le.nysite) then
    iysucc=iysite+1
 else
    iysucc=1
 endif
 RETURN
 END FUNCTION iysucc

 END MODULE site_finder_module 
