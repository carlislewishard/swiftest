!**********************************************************************************************************************************
!
!  Unit Name   : symba_regime
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Determine the collisional regime of two colliding bodies.
!
!  Input
!    Arguments : 
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : 
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_regime(Mcenter, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, mtiny)
!
!  Notes       : Current version requires all values to be converted to SI units prior to calling the function
!                 Reference:
!                 Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
!                    Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
!
!
!**********************************************************************************************************************************
SUBROUTINE symba_regime(Mcenter, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, mtiny)

! Modules
     USE swiftest
     USE module_symba
     USE module_helio
     USE module_nrutil
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => symba_regime
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT)                 :: regime
     REAL(DP), INTENT(OUT)                     :: Mlr, Mslr
     REAL(DP), INTENT(IN)                      :: Mcenter, m1, m2, rad1, rad2, den1, den2, mtiny 
     REAL(DP), DIMENSION(:), INTENT(IN)     :: xh1, xh2, vb1, vb2

! Internals
     REAL(DP)                      :: a1, alpha, Aint, b, bcrit, E, fgamma, l, Lint, mu, phi, theta
     REAL(DP)                      :: QR, QRD_pstar, QR_erosion, QR_supercat
     REAL(DP)                      :: vcr, verosion, vescp, vhill, vimp, vsupercat
     REAL(DP)                      :: mint, mtot
     REAL(DP)                      :: Rp, Rhill 
! Constants
     INTEGER(I4B)                  :: N1 = 1  !number of objects with mass equal to the largest remnant from LS12
     INTEGER(I4B)                  :: N2 = 2  !number of objects with mass larger than second largest remnant from LS12
     !INTEGER(I4B)                  :: N1g = 2  !number of objects with mass equal to the largest remnant from LS12 if Mp = Mtarg
     !INTEGER(I4B)                  :: N2g = 4  !number of objects with mass larger than second largest remnant from LS12 if Mp = Mtarg
     REAL(DP)                      :: density1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
     REAL(DP)                      :: c_star = 1.8_DP !3.0 #3.0# #5#1.8 #1.8 #Measure of dissipation of energy within the target (Chambers frag.f90)
     REAL(DP)                      :: mu_bar = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
     REAL(DP)                      :: beta = 2.85_DP !slope of SFD for remnants from LS12 2.85
     REAL(DP)                      :: c1 = 2.43_DP !Ls12 constants
     REAL(DP)                      :: c2 = -0.0408_DP !Ls12 constants
     REAL(DP)                      :: c3 = 1.86_DP !Ls12 constants
     REAL(DP)                      :: c4 = 1.08_DP !Ls12 constants
     REAL(DP)                      :: c5 = 2.5_DP !Ls12 constants
     REAL(DP)                      :: crufu = (2.0_DP-3.0_DP*0.36_DP) ! central potential variable from Rufu et al. 2019
     real(DP)                      :: mresidual

! Executable code
      vimp = norm2(vb2(:) - vb1(:))
      b = calc_b(xh2, vb2, xh1, vb1)
      l = (rad1 + rad2) * (1 - b)
      E = (norm2(vb1)**2) / 2.0_DP - GC * Mcenter / norm2(xh1)
      a1 = - GC * Mcenter / 2.0_DP / E
      mtot = m1 + m2 
      mu = (m1 * m2) / mtot
      IF (l < 2 * rad2) THEN
            !Calculate mint
            Phi = 2 * acos((l - rad2) / rad2)
            Aint = (rad2**2) * (PI - (Phi - sin(Phi)) / 2.0_DP)
            Lint = 2 * sqrt(rad2**2 - (rad2 - l / 2.0_DP) ** 2) 
            mint = Aint * Lint  ![kg]
            alpha = (l**2) * (3 * rad2 - l) / (4 * (rad2**3))
      ELSE
           alpha = 1.0_DP
           mint = m2
      END IF 
      Rp = (3 * (m1 / den1 + alpha * m2 / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! (Mustill et al. 2019)
     !Calculate vescp
      vescp = sqrt(2 * GC * (mtot) / (Rp)) !Mustill et al. 2018 Eq 6 
     !Calculate Rhill
      Rhill = a1 * (m1 / 3.0_DP / (Mcenter + m1))**(1.0_DP/3.0_DP)
     !Calculate Vhill
      if ((rad2 + rad1) < Rhill) then 
        vhill = sqrt(2 * GC * m1 * ((Rhill**2 - Rhill * (rad1 + rad2)) / &
          (Rhill**2 - 0.5_DP * (rad1 + rad2)**2)) / (rad1 + rad2))
      else
        vhill = vescp
      end if 
     !Calculate QR_pstar
      QRD_pstar = calc_QRD_pstar(m1, m2, alpha) * (vhill / vescp)**crufu !rufu et al. eq (3)
     !Calculate verosion
      QR_erosion = 2 * (1.0_DP - m1 / mtot) * QRD_pstar
      verosion = (2* QR_erosion * mtot / mu)** (1.0_DP / 2.0_DP)
      QR = mu*(vimp**2) / mtot / 2.0_DP
      !QRD_lr = calc_QRD_lr(m2, m1, mint)
     !Calculate Mass largest remnant Mlr 
      Mlr = (1.0_DP - QR / QRD_pstar / 2.0_DP) * mtot  ! [kg] #(Eq 5)
     !Calculate vsupercat
      QR_supercat = 1.8_DP * QRD_pstar
      vsupercat = sqrt(2 * QR_supercat * mtot / mu)
     !Calculate Vcr
      fgamma = (m1 - m2) / mtot
      theta = 1.0_DP - b
      vcr = vescp * (c1 * fgamma * theta**c5 + c2 * fgamma + c3 * theta**c5 + c4)
      bcrit = rad1 / (rad1 + rad2)

      IF ((m1 < 1*MTINY).OR.(m2 < 1*MTINY)) THEN 
        regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
        Mlr = mtot
        Mslr = 0.0_DP
        WRITE(*,*) "FORCE MERGE"
      

      ELSE 

        IF( vimp < vescp) THEN
          regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
          Mlr = mtot
          Mslr = 0.0_DP
        ELSE IF (vimp < verosion) THEN 
          IF (b < bcrit) THEN
            regime = COLLRESOLVE_REGIME_MERGE !partial accretion regime"
            Mlr = mtot
            Mslr = 0.0_DP
          ELSE IF ((b > bcrit) .AND. (vimp < vcr)) THEN
            regime = COLLRESOLVE_REGIME_MERGE ! graze and merge
            Mlr = mtot
            Mslr = 0.0_DP
          ELSE
            Mlr = m1
            Mslr = calc_QRD_rev(m2,m1,mint,den1,den2,vimp)
            regime = COLLRESOLVE_REGIME_HIT_AND_RUN !hit and run
          END IF 
        ELSE IF (vimp > verosion .AND. vimp < vsupercat) THEN
          IF ((m2 < 0.001_DP * m1)) THEN 
            regime = COLLRESOLVE_REGIME_MERGE !cratering regime"
            Mlr = mtot
            Mslr = 0.0_DP
          ELSE 
            Mslr = (mtot * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / mtot)))) / (N2 * beta)  ! (Eq 37)
            regime = COLLRESOLVE_REGIME_DISRUPTION !disruption
          END IF 
        ELSE IF (vimp > vsupercat) THEN 
          Mlr = mtot * (0.1_DP * ((QR / (QRD_pstar * 1.8_DP))**(-1.5_DP)))     !Eq (44)
          Mslr = mtot * (3.0_DP - beta) * (1.0_DP - N1 * Mlr / mtot) / (N2 * beta)  ! (Eq 37)
          regime = COLLRESOLVE_REGIME_SUPERCATASTROPHIC ! supercatastrophic
        ELSE 
          WRITE(*,*) "Error no regime found in symba_regime"
        END IF 
      END IF 

      mresidual = mtot - Mlr - Mslr
      if (mresidual < 0.0_DP) then ! Prevents final masses from going negative
         Mlr = Mlr + mresidual
      end if
         
    RETURN 


! Functions
contains
function calc_QRD_pstar(Mtarg,Mp,alpha) result(ans)
   implicit none
   real(DP),intent(in) :: Mtarg, Mp, alpha
   real(DP)            :: QRD_star1, mu_alpha, mu, QRD_star, QRD_pstar
   real(DP)            :: ans
   ! calc mu, mu_alpha
   mu = (Mtarg * Mp) / (Mtarg + Mp)  ! [kg]
   mu_alpha = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)  ! [kg]
   ! calc QRD_star1
   QRD_star1 = (c_star * 4 * PI * density1 * GC * Rp**2) / 5.0_DP
   ! calc QRD_star
   QRD_star = QRD_star1 * (((Mp / Mtarg + 1.0_DP)**2) / (4 * Mp / Mtarg))**(2.0_DP / (3.0_DP * mu_bar) - 1.0_DP)  !(Eq 23)
   ! calc QRD_pstar, V_pstar
   QRD_pstar = ((mu / mu_alpha)**(2.0_DP - 3.0_DP * mu_bar / 2.0_DP)) * QRD_star  ! (Eq 15)
   
   ans = QRD_pstar
   return
end function calc_QRD_pstar

function calc_QRD_rev(Mp,Mtarg,mint,den1,den2, vimp) result(ans)
   implicit none
   real(DP),intent(in) :: Mp, Mtarg, mint, den1, den2, vimp
   real(DP) :: ans, Mtot_rev, mu_rev, gamma_rev, QRD_star1, QRD_star, mu_alpha_rev
   real(DP) :: QRD_pstar, RC1, QR_rev, QRD_pstar_rev, Mslr, QR_supercat_rev
   ! calc Mtlr, RC1, mu, gammalr
   Mtot_rev =  mint + Mp
   RC1 = (3 * (mint / den1 + Mp / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! [m] Mustill et al 2018
   mu_rev = (mint * Mp) / Mtot_rev ! [kg] Eq 49 LS12
   mu_alpha_rev = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)
   gamma_rev = mint / Mp ! Eq 50 LS12
   !calc QR_rev
   QR_rev = mu_rev * (vimp**2) / (2 * Mtot_rev)
   ! calc QRD_star1, V_star1
   QRD_star1 = (c_star * 4 * PI * Mtot_rev * GC ) / RC1 / 5.0_DP
   ! calc QRD_pstar_rev
   QRD_star = QRD_star1 * (((gamma_rev + 1.0_DP)**2) / (4 * gamma_rev)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP) !(Eq 52)
   QRD_pstar = QRD_star * ((mu_rev / mu_alpha_rev)**(2.0_DP - 3.0_DP * mu_bar / 2.0_DP))
   QRD_pstar_rev = QRD_pstar * (vhill / vescp)**crufu !rufu et al. eq (3)
   !calc QR_supercat_rev
   QR_supercat_rev = 1.8_DP * QRD_pstar_rev 
   !V_supercat_rev = ( 2.0_DP * QR_supercat_rev * Mtot_rev / mu_rev ) ** (1.0_DP / 2.0_DP)
   if (QR_rev > QR_supercat_rev ) then 
      Mslr = Mtot_rev * (0.1_DP * ((QR_rev / (QRD_pstar_rev * 1.8_DP))**(-1.5_DP)))     !Eq (44)
   else if ( QR_rev < QRD_pstar_rev ) then 
      Mslr = mp 
   else 
      Mslr = (1.0_DP - QR_rev / QRD_pstar_rev / 2.0_DP) * (Mtot_rev)  ! [kg] #(Eq 5)
   end if 

   if ( Mslr > mp ) Mslr = mp !Check conservation of mass
   ans = Mslr

   return
end function calc_QRD_rev

function calc_b(proj_pos, proj_vel, targ_pos, targ_vel) result(ans)
  implicit none
  real(DP), intent(in), dimension(:) :: proj_pos, proj_vel, targ_pos, targ_vel
  real(DP)                           :: angle, ans 
  real(DP), dimension(NDIM)          :: imp_vel, distance         

    imp_vel = proj_vel - targ_vel
    distance = proj_pos - targ_pos
    angle = acos(dot_product(imp_vel, distance) / norm2(imp_vel) / norm2(distance))      
    ans = sin(angle)
  return 
end function calc_b

! function calc_b(Mp_pos, Mp_vel, Mp_r, Mtarg_pos, Mtarg_vel, Mtarg_r) result(b)
!    implicit none
!    real(DP), intent(in), DIMENSION(3) :: Mp_pos, Mp_vel, Mtarg_pos, Mtarg_vel
!    real(DP), intent(in) :: Mp_r, Mtarg_r
!    real(DP) :: h_sq, b, dvel_sq
!    real(DP), DIMENSION(3) :: dpos, dvel, h

!    dpos(1) = mtarg_pos(1) - mp_pos(1)
!    dpos(2) = mtarg_pos(2) - mp_pos(2)
!    dpos(3) = mtarg_pos(3) - mp_pos(3)

!    dvel(1) = mtarg_vel(1) - mp_vel(1)
!    dvel(2) = mtarg_vel(2) - mp_vel(2)
!    dvel(3) = mtarg_vel(3) - mp_vel(3)

!    h(1) = (dpos(2) * dvel(3)) - (dpos(3) * dvel(2))
!    h(2) = (dpos(3) * dvel(1)) - (dpos(1) * dvel(3))
!    h(3) = (dpos(1) * dvel(2)) - (dpos(2) * dvel(1))

!    h_sq = (h(1) * h(1)) + (h(2) * h(2)) + (h(3) * h(3))
!    dvel_sq = (dvel(1) * dvel(1)) + (dvel(2) * dvel(2)) + (dvel(3) * dvel(3))

!    b = (h_sq / (((Mp_r + Mtarg_r) ** 2.0_DP) * dvel_sq)) ** (1.0_DP / 2.0_DP)
!    print(b,"b")
!    return
! end function calc_b

END SUBROUTINE symba_regime
!**********************************************************************************************************************************
!
!  Author(s)   : C.Wishard and J.Pouplin
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
