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


 MODULE cubic_eos_module 
 USE variable
 IMPLICIT NONE

 CONTAINS      


 SUBROUTINE calc_mix_terms(mixing_rule,pure_term,y,int_param,mix_term)      !This subroutine calculate the mix terms.
 CHARACTER*2, INTENT(IN)       :: mixing_rule 
 DOUBLE PRECISION, INTENT(IN)  :: pure_term(ncomp),y(ncomp)       !y is gas phase mole fraction
 DOUBLE PRECISION, INTENT(IN)  :: int_param(1:ncomp,1:ncomp)      !Binary interaction parameter
 DOUBLE PRECISION, INTENT(OUT) :: mix_term(1:ncomp,1:ncomp)

 INTEGER                       :: icomp, jcomp                    !dummy counter variables

 Do icomp=1,ncomp
    mix_term(icomp,icomp)=pure_term(icomp)
    Do jcomp=icomp+1,ncomp                                        ! so icomp < jcomp Only part of the matrix is calculated
       select case (mixing_rule)
       case("AM")
                 mix_term(icomp,jcomp)=0.5d0*(pure_term(icomp)+pure_term(jcomp))*(1-int_param(icomp,jcomp))         
       case("GM")
                 mix_term(icomp,jcomp)=sqrt(pure_term(icomp)*pure_term(jcomp))*(1-int_param(icomp,jcomp))
       end select
       mix_term(jcomp,icomp)=mix_term(icomp,jcomp)               !By symmetry property of mix term
    Enddo 
 Enddo
 RETURN
 END SUBROUTINE calc_mix_terms


 SUBROUTINE calc_tot_mix_and_cross_term(mix_term,y,tot_cross,tot_mix)  ! This subroutine calculates total mix term and cross term for each components
 DOUBLE PRECISION, INTENT(IN) :: mix_term(1:ncomp,1:ncomp),y(ncomp)
 DOUBLE PRECISION, INTENT(OUT) :: tot_cross(ncomp),tot_mix
 INTEGER          :: icomp, jcomp
 
 tot_mix=0.0d0                                                     !Initialization of variables
 Do icomp=1,ncomp
    tot_cross(icomp)=0.0d0
 Enddo

 Do icomp=1,ncomp                                                  !calculation of sums
    Do jcomp=1,ncomp
       tot_cross(icomp)=tot_cross(icomp)+y(jcomp)*mix_term(icomp,jcomp)
       tot_mix=tot_mix+y(icomp)*y(jcomp)*mix_term(icomp,jcomp)
    Enddo
 Enddo
 RETURN
 END SUBROUTINE calc_tot_mix_and_cross_term



 SUBROUTINE eos_2param_solve(eos_type,across,bcross,amix,bmix,phi) !Subroutine to calculate Redlick-Kwong equation type fugacities
 CHARACTER*3,INTENT(IN)        :: eos_type
 DOUBLE PRECISION, INTENT(IN)  :: across(ncomp),bcross(ncomp)      !Input arrays are critical pressure and temperature for each component and composition
 DOUBLE PRECISION, INTENT(IN)  :: amix,bmix
 DOUBLE PRECISION, INTENT(OUT) :: phi(ncomp)                       !Outputs are the fugacity coefficient
 DOUBLE PRECISION              :: lnphi                            ! temporary variable to store lnphi 
 DOUBLE PRECISION              :: q,r,s,z                          !Coefficients of cubic equation Z^3+q*Z^2+r*Z+s==0; Z=z=Compression factor
 DOUBLE PRECISION              :: prefactor,rt2,log_arg_pr         !prefactor of term1, sqrt_of_2, argument for log function of PR eos
 DOUBLE PRECISION :: term1(ncomp),term2(ncomp),term3(ncomp)        !component dependent terms
 INTEGER :: ic                                                     !Counter 

 select case(eos_type)

 case("VdW")
 q = -1.0*(bmix+1)
 r = amix
 s = -1.0*amix*bmix
 z = cube_root(q,r,s)
! write(*,*) "Z=", z
 do ic=1,ncomp
    term1(ic) = (amix-2.0d0*across(ic))/z
    term2(ic) = 2.0d0*(bcross(ic)-bmix)/(z-bmix)
    term3(ic) = 0.0d0
 enddo

 case("RKX")
 q = -1.0d0
 r =  amix-bmix-bmix*bmix
 s = -1.0d0*amix*bmix  
 z =  cube_root(q,r,s)
! write(*,*) "Z=", z      
 do ic=1,ncomp
    prefactor=(2.0d0*amix*bcross(ic)-2.0d0*across(ic)*bmix-amix*bmix)/(bmix*bmix)
    term1(ic) =  prefactor*dlog(1+bmix/z)
    term2(ic) =  2.0d0*(bcross(ic)-bmix)/(z-bmix)
    term3(ic) = -2.0d0*(bcross(ic)-bmix)*amix/(bmix*(z+bmix)) 
!*   write(*,*)"(term2+term3)/(Z-1)=", (term2(ic)+term3(ic))/(z-1) !To check the formula
!*   write(*,*) "prefactor=",prefactor
 enddo

 case("PRX")
 rt2=sqrt(2.0d0)
 q = bmix-1
 r = amix-3.0d0*bmix*bmix-2.0*bmix
 s = bmix**3+bmix**2-amix*bmix
 z = cube_root(q,r,s)
! write(*,*) "Z=", z
! write(*,*) "amix=",amix,"bmix=",bmix
 do ic=1,ncomp
    prefactor = (amix*bcross(ic)-across(ic)*bmix-0.5d0*amix*bmix)/(rt2*bmix*bmix)
    log_arg_pr = (z+bmix+rt2*bmix)/(z+bmix-rt2*bmix)
    term1(ic) = prefactor*dlog(log_arg_pr)
    term2(ic) = 2.0d0*(bcross(ic)-bmix)/(z-bmix)
    term3(ic) = 2.0d0*z*amix*(bmix-bcross(ic))/(bmix*(z*z+2.0d0*z*bmix-bmix*bmix))
!*   write(*,*) "(term2+term3)/(z-1)=", (term2(ic)+term3(ic))/(z-1)             !To check the formula
!*   write(*,*) 'prefactor=',prefactor,'across=',across(ic),'bcross',bcross(ic) !To check formula
 enddo
 end select

 do ic=1,ncomp
    lnphi = z-1.0d0-dlog(z-bmix)+term1(ic)+term2(ic)+term3(ic)
    phi(ic) = dexp(lnphi)
 enddo
 RETURN
 END SUBROUTINE eos_2param_solve                                   !End of subroutine eos_2param_solve

     
 SUBROUTINE calculate_fugacity(p,t,y,phi)                     !Subrutine that calculates fugacity coefficients for each component according to given EOS
 DOUBLE PRECISION, INTENT(IN)  :: p,t                             !pressure, temperature
 DOUBLE PRECISION, INTENT(IN)  :: y(ncomp)                        !gas composition
 DOUBLE PRECISION, INTENT(OUT) :: phi(ncomp)                      ! Fugacity coefficients
 DOUBLE PRECISION :: pr(ncomp),tr(ncomp),acc_fac(ncomp)           !reduced pressure, reduced temperature, accentric_factor 
 DOUBLE PRECISION :: aa(ncomp),bb(ncomp),cc(ncomp)                ! VdW a and b parameter and c parameter for PTV for each components
 DOUBLE PRECISION :: a(1:ncomp,1:ncomp), b(1:ncomp,1:ncomp)       !VdW a and b parameters for each component as well as the mixed terms
 DOUBLE PRECISION :: c(1:ncomp,1:ncomp)                           !PTV c terms for each component as well as mixed terms
 DOUBLE PRECISION :: across(ncomp),bcross(ncomp),ccross(ncomp)    !total cross terms for each component
 DOUBLE PRECISION :: amix,bmix,cmix                               !total mixed terms 
 DOUBLE PRECISION :: temp_func, m, k0, k1                         !temporary variable to stor temp_function related numbers
 DOUBLE PRECISION :: prg1,prg2,prg_tr_exp                         !PRG related some temporary variables
 INTEGER :: icomp

 do icomp=1,ncomp
    pr(icomp)=p/comp_information(icomp)%pc
    tr(icomp)=t/comp_information(icomp)%tc
    acc_fac(icomp)=comp_information(icomp)%acc_fac
 enddo

 select case(eos)
 
 case("Ideal")
 write(*,'(A)') "Calculating fugacity coefficients according to Ideal EOS."
 do icomp=1,ncomp
    phi(icomp)=1.0d0
 enddo

 case("VdW") 
 write(*,'(A)') "Calculating fugacity coefficients according to VdW EOS."
 do icomp=1,ncomp
    aa(icomp)=(27.0d0*pr(icomp))/(64.0d0*tr(icomp)**2)
    bb(icomp)=pr(icomp)/(8.0d0*tr(icomp))
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("VdW",across,bcross,amix,bmix,phi)
 
 case("RK")
 write(*,'(A)') "Calculating fugacity coefficients according to RK EOS."
  do icomp=1,ncomp
    aa(icomp)=(0.42747d0*pr(icomp))/(tr(icomp)**2.5)
    bb(icomp)=0.08664d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("RKX",across,bcross,amix,bmix,phi)     

 case("SRK")
 write(*,'(A)') "Calculating fugacity coefficients according to SRK EOS."
 do icomp=1,ncomp
    m = 0.480d0+1.574d0*acc_fac(icomp)-0.176d0*acc_fac(icomp)**2
    temp_func = (1+m*(1.0d0-dsqrt(tr(icomp))))**2
    aa(icomp)=(0.42747d0*pr(icomp)*temp_func)/(tr(icomp)**2)
    bb(icomp)=0.08664d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("RKX",across,bcross,amix,bmix,phi)

!**************FOR DEBUG*******************************************
!*  write(*,*) "pure_prefactor:", amix/bmix
!*  do icomp=1,ncomp
!*     write(*,*)"(b/bmix-1)", bb(icomp)/bmix-1
!*     write(*,*) "additional prefactor:",2.0d0*dsqrt(aa(icomp)/amix)-bb(icomp)/bmix
!*  enddo
!*******************************************************************

 case("PR")
 write(*,'(A)') "Calculating fugacity coefficients according to PR EOS"
 do icomp=1,ncomp
    m=0.37464d0+1.54226*acc_fac(icomp)-0.26992*acc_fac(icomp)**2
    temp_func = (1+m*(1.0d0-dsqrt(tr(icomp))))**2
    aa(icomp)=(0.45724d0*pr(icomp)*temp_func)/(tr(icomp)**2)
    bb(icomp)=0.07780d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("PRX",across,bcross,amix,bmix,phi)

!**************FOR DEBUG*******************************************
!* write(*,*) "pure_prefactor=",amix/(2.0d0*dsqrt(2.0d0)*bmix)
!* do icomp=1,ncomp
!*    write(*,*)"(b/bmix-1)", bb(icomp)/bmix-1
!*    write(*,*) "additional_prefactor",2.0d0*across(icomp)/amix-bb(icomp)/bmix
!* enddo
!*******************************************************************


 case("PR78")
 write(*,'(A)') "Calculating fugacity coefficients according to PR78 EOS"
 do icomp=1,ncomp
    if(acc_fac(icomp).le.0.491) then
      m=0.37464d0+1.54226*acc_fac(icomp)-0.26992*acc_fac(icomp)**2
    else
      m=0.379642+1.48503*acc_fac(icomp)-0.164423*acc_fac(icomp)**2+0.016666*acc_fac(icomp)**3
    endif
    temp_func = (1+m*(1.0d0-dsqrt(tr(icomp))))**2                 
    aa(icomp)=(0.45724d0*pr(icomp)*temp_func)/(tr(icomp)**2)
    bb(icomp)=0.07780d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("PRX",across,bcross,amix,bmix,phi)

 case("PR80")
 write(*,'(A)') "Calculating fugacity coefficients according to PR80 EOS"
 do icomp=1,ncomp
    temp_func = (1.0085677d0+0.82154*(1.0d0-dsqrt(tr(icomp))))**2  
    aa(icomp)=(0.45724d0*pr(icomp)*temp_func)/(tr(icomp)**2)
    bb(icomp)=0.07780d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("PRX",across,bcross,amix,bmix,phi)

 case("PRG")
 write(*,'(A)') "Calculating fugacity coefficients according to PRG EOS"
 do icomp=1,ncomp
    prg1=2.0d0+0.836d0*tr(icomp)
    prg_tr_exp=0.134+0.508*acc_fac(icomp)-0.0467*acc_fac(icomp)**2
    prg2=1.0d0-tr(icomp)**prg_tr_exp
    temp_func = exp(prg1*prg2)
    aa(icomp)=(0.45724d0*pr(icomp)*temp_func)/(tr(icomp)**2)
    bb(icomp)=0.07780d0*pr(icomp)/tr(icomp)
 enddo
 call calc_mix_terms("GM",aa,y,int_param_a,a)
 call calc_mix_terms("AM",bb,y,int_param_b,b)
 call calc_tot_mix_and_cross_term(a,y,across,amix)
 call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
 call eos_2param_solve("PRX",across,bcross,amix,bmix,phi)

! case("gPRSV")
! write(*,'(A)') "Calculating fugacity coefficients according to gPRSV EOS"
! do icomp=1,ncomp
!    k0=0.379368d0+1.459994*zc(icomp)-0.125563*acc_fac(icomp)**2
!    k1=0.599529d0-1.952083*zc(icomp)+0.080764*acc_fac(icomp)-0.209272*acc_fac(icomp)**2
!    m=k0+k1*(1+dsqrt(tr(icomp)))*(0.7d0-tr(icomp))
!    temp_func = (1.0d0+m*(1.0d0-dsqrt(tr(icomp))))**2
!    aa(icomp)=(0.45724d0*pr(icomp)*temp_func)/(tr(icomp)**2)
!    bb(icomp)=0.07780d0*pr(icomp)/tr(icomp)
! enddo
! call calc_mix_terms("GM",aa,y,int_param_a,a)
! call calc_mix_terms("AM",bb,y,int_param_b,b)
! call calc_tot_mix_and_cross_term(a,y,across,amix)
! call calc_tot_mix_and_cross_term(b,y,bcross,bmix)
! call eos_2param_solve("PRX",across,bcross,amix,bmix,phi)
 
 end select
 RETURN
 END SUBROUTINE calculate_fugacity

       
 FUNCTION cube_root(a,b,c)                                        !See, Numerical recipes in Fortran, section 5.6, for the algorithm.
 IMPLICIT NONE
 DOUBLE PRECISION :: cube_root
 DOUBLE PRECISION, PARAMETER :: PI=3.1415926535897932
 DOUBLE PRECISION, PARAMETER :: EPS = 1.0E-17
 DOUBLE PRECISION, INTENT(IN) :: a,b,c
 DOUBLE PRECISION :: q,qcube,qsqrt,r,rsq,aa,bb,x1,x2,x3
 DOUBLE PRECISION :: theta,cube_root_old,root_diff
 DOUBLE PRECISION :: cube_func_val, cube_func_deriv_val
 INTEGER          :: newton_iteration_step
 REAL             :: newton_start_time, newton_end_time

 root_diff=1.0d0
 q=(a*a-3.0d0*b)/9.0d0
 r=(2.0d0*a*a*a-9.0d0*a*b+27.0d0*c)/54.0d0
 qcube=q*q*q
 rsq=r*r

 if(rsq.lt.qcube) then
!****   write(*,*) "R2 is less thant Q3."
   theta=dacos(r/(q**1.5))
   qsqrt=dsqrt(q)
   x1=-2.0d0*qsqrt*dcos(theta/3.0d0)-a/3.0d0
   x2=-2.0d0*qsqrt*dcos((theta+2.0d0*pi)/3.0d0)-a/3.0d0
   x3=-2.0d0*qsqrt*dcos((theta-2.0d0*pi)/3.0d0)-a/3.0d0

   cube_root=max(x1,x2,x3)
 else
!****   write(*,*) "R2 is greater than or equalto Q3."
   aa=-sign(1.0d0,r)*(abs(r)+sqrt(rsq-qcube))**(1.0d0/3.0d0)
   if(dabs(aa).eq.0.0d0) then
     bb=0.0d0
   else
     bb=q/aa
   endif
   
   cube_root=aa+bb-a/3.0d0
 endif
 
 newton_iteration_step=0
 call cpu_time(newton_start_time)
 do while ((root_diff.gt.EPS).and.(newton_iteration_step.lt.1000000000))                    !Additional Newton step to refine cuberoot
          cube_func_val=cube_root*cube_root*cube_root+a*cube_root*cube_root+b*cube_root+c   !original cubic function whose roots to be found
          cube_func_deriv_val=3.0d0*cube_root*cube_root+2.0d0*a*cube_root+b                 ! derivative of the cubic function
          cube_root_old=cube_root
          cube_root=cube_root_old-cube_func_val/cube_func_deriv_val                         !Newton iteration step

          root_diff=Abs(cube_root-cube_root_old)
          newton_iteration_step=newton_iteration_step+1
!****        write(*,*) "Newton iteration no",newton_iteration_step   
 enddo   
!**** write(*,*) "Total No of Newton iteration =", newton_iteration_step
 call cpu_time(newton_end_time)
!**** write(*,*) "Newton runtime(s)=", newton_end_time-newton_start_time
 RETURN
 END FUNCTION cube_root

 END MODULE cubic_eos_module

