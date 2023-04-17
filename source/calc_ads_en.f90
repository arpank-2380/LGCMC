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


 MODULE calc_ads_en
 USE variable
 USE site_finder_module, ONLY : site_finder

 CONTAINS
 FUNCTION  ads_en(ixsite,iysite,isite_type,icomp)
 DOUBLE PRECISION :: ads_en
 INTEGER, INTENT(IN) :: ixsite,iysite,isite_type,icomp
! LOGICAL :: constraint_satisfied 
 
 if (icomp.eq.0) then
    ads_en=0.0d0
    RETURN

 else if (.not.constraint_present(icomp,isite_type)) then
    ads_en=h(icomp,isite_type)
    RETURN
 end if


 if(constraint_satisfied(ixsite,iysite,isite_type,icomp)) then
   ads_en=h(icomp,isite_type)
 else
   ads_en=h_constrain(icomp,isite_type)
 endif

!****************************DEBUG********************************************
!*******      if(ads_en.ge.h_cons_min) then
!******        Print *, 'Constraints FALSE',ixsite,iysite,isite_type,
!******     x  icomp,confg(ixsite,iysite,1)
!******      endif
!*****************************************************************************
 RETURN
 END FUNCTION ads_en


 FUNCTION constraint_satisfied(ixsite,iysite,isite_type,icomp)
 INTEGER, INTENT(IN) :: ixsite,iysite,isite_type,icomp
 INTEGER :: ix,iy                                  !direction of constraint
 INTEGER :: jxsite,jysite,jsite_type,jcomp,kcomp      

 LOGICAL :: constraint_satisfied
 LOGICAL :: constraint_satisfied_temp

! Declaration of counter variables
 INTEGER icons,ii,jj

  constraint_satisfied=.TRUE.

 if(icomp.eq.0) then
   RETURN
 endif

 if(.not.constraint_present(icomp,isite_type)) then
    RETURN
 endif
 
 do icons=1,nconstraints(icomp,isite_type)
    jsite_type=cons(icomp,isite_type,icons)%site
    ix=cons(icomp,isite_type,icons)%x
    iy=cons(icomp,isite_type,icons)%y
     
    call site_finder(ixsite,iysite,ix,iy,jxsite,jysite)
     
    jcomp=confg(jxsite,jysite,jsite_type)

    constraint_satisfied_temp=.FALSE.     
! if switch is true then it will try to search the allowed component
      if(jcomp.ne.0) then 
         if(cons(icomp,isite_type,icons)%comp_list(jcomp).eqv.cons(icomp,isite_type,icons)%switch) then
            constraint_satisfied_temp=.TRUE.
!****      print *,'True encountered'
         endif
          
      else
         if(.not.cons(icomp,isite_type,icons)%switch) then
             constraint_satisfied_temp=.TRUE.
!****      Print *, 'True encountered'
         endif
      endif 

    constraint_satisfied=constraint_satisfied.and.constraint_satisfied_temp

    if(.not.constraint_satisfied) then
!****         Print *, "Constraint become False"
!****         Print *,ixsite,iysite,isite_type,icomp,confg(ixsite,iysite,1)
      exit
    endif
 enddo

!****      if (constraint_satisfied) then
!****      print *, 'constraint satisfied?',constraint_satisfied
!****      endif
   
 RETURN
 END FUNCTION constraint_satisfied
 
 END MODULE calc_ads_en 
