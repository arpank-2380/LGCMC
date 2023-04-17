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


 MODULE calc_lat_int_en
 USE site_finder_module, ONLY : site_finder
 IMPLICIT NONE
 INTEGER, PRIVATE :: iixsite, iiysite,iisite_type,iicomp          !Dummy variables are changed into accessible one
 INTEGER, PRIVATE :: jxsite,jysite,jsite_type,jcomp
 INTEGER, PRIVATE :: ix,iy

 CONTAINS

 FUNCTION  lat_int_en(ixsite,iysite,isite_type,icomp)
 USE variable, ONLY : xcut,ycut,confg,nsite_type,read_lat_int
 DOUBLE PRECISION lat_int_en,jpair
 INTEGER, INTENT(IN) :: ixsite,iysite,isite_type,icomp
 INTEGER xmin,xmax,ymin,ymax                           !xcut ycut range

 lat_int_en = 0.0d0                       !Initialization of lat_int_en

 if (icomp.eq.0) then                     !If vacuum then there would be no lateral interactions
    RETURN                                !i.e it will return the value 0 then
 endif

 xmin = -1*xcut
 ymin = -1*ycut
 xmax = xcut
 ymax = ycut

!Assignment of dummy variable into global variables
 iixsite=ixsite
 iiysite=iysite
 iisite_type=isite_type
 iicomp=icomp

 Do ix=xmin,xmax
    Do iy=ymin,ymax
       call site_finder(ixsite,iysite,ix,iy,jxsite,jysite)  
       Do jsite_type=1,nsite_type
          jcomp=confg(jxsite,jysite,jsite_type)  
          if(jcomp.eq.0) then
            cycle
          
          else if(read_lat_int(icomp,isite_type,ix,iy,jcomp,jsite_type)) then 
               call pair_int_en(jpair)
              lat_int_en= lat_int_en+jpair
!****      write(*,'(6I3,2x,G10.5)') ix,iy,jsite_type,jxsite,jysite,jcomp,jpair
          endif       
       Enddo
    Enddo
 Enddo

 RETURN
 END FUNCTION lat_int_en

 SUBROUTINE pair_int_en(jpair) 
 USE variable, ONLY : asym,j
 DOUBLE PRECISION :: jpair,j_pair_asym

 if(asym(iicomp,iisite_type,ix,iy,jcomp,jsite_type)) then
    call lattice_assymmetry(j_pair_asym)
    jpair=j_pair_asym
 else
    jpair=j(iicomp,iisite_type,ix,iy,jcomp,jsite_type)
 endif   
 RETURN
 END SUBROUTINE pair_int_en

 SUBROUTINE lattice_assymmetry(j_pair_asym)
 USE variable, ONLY : nxsite, nysite, j, jb
 DOUBLE PRECISION :: j_pair_asym
 LOGICAL :: sublattice_odd
 INTEGER :: shift_frm_reference
 
 INTEGER :: numerator,denominator,increment
 
 if(Abs(ix).ge.Abs(iy)) then
   numerator=iixsite
   denominator=Abs(ix)
   increment=nxsite
 else
   numerator=iiysite
   denominator=Abs(iy)
   increment=nysite
 endif
 
 if (denominator.eq.1) then
    shift_frm_reference=numerator
 else  
    Do while (Mod(numerator,denominator).ne.1) 
       numerator=numerator+increment
!****      write(*,*) numerator
    Enddo
       shift_frm_reference=numerator/denominator
 endif
      
 if(Mod(shift_frm_reference,2).eq.1) then
   j_pair_asym=j(iicomp,iisite_type,ix,iy,jcomp,jsite_type)
 else
   j_pair_asym=jb(iicomp,iisite_type,ix,iy,jcomp,jsite_type)
 endif

 RETURN
 END SUBROUTINE lattice_assymmetry

 END MODULE calc_lat_int_en 
