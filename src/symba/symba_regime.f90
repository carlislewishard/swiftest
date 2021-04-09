subroutine symba_regime(Mcb, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, mlr, mslr, mtiny)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Determine the collisional regime of two colliding bodies. 
   !! Current version requires all values to be converted to SI units prior to calling the function
   !!       Reference:
   !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
   !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
   !!       Mustill, A.J., Davies, M.B., Johansen, A., 2018. The dynamical evolution of transiting planetary systems including 
   !!          a realistic collision prescription. Mon Not R Astron Soc 478, 2896–2908. https://doi.org/10.1093/mnras/sty1273
   !!       Rufu, R., Aharonson, O., 2019. Impact Dynamics of Moons Within a Planetary Potential. J. Geophys. Res. Planets 124, 
   !!          1008–1019. https://doi.org/10.1029/2018JE005798
   !!
   use swiftest
   use module_symba
   use module_helio
   use module_nrutil
   use module_swiftestalloc
   use module_interfaces, EXCEPT_THIS_ONE => symba_regime
   implicit none

! Arguments
   integer(I4B), intent(out)         :: regime
   real(DP), intent(out)          :: mlr, mslr
   real(DP), intent(in)           :: Mcb, m1, m2, rad1, rad2, den1, den2, mtiny 
   real(DP), dimension(:), intent(in)   :: xh1, xh2, vb1, vb2
! Constants
   integer(I4B), parameter :: N1 = 1  !number of objects with mass equal to the largest remnant from LS12
   integer(I4B), parameter :: N2 = 2  !number of objects with mass larger than second largest remnant from LS12
   real(DP), parameter   :: DENSITY1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
   real(DP), parameter   :: C_STAR = 1.8_DP !3.0 #3.0# #5#1.8 #1.8 #measure of dissipation of energy within the target (chambers frag.f90)
   real(DP), parameter   :: MU_BAR = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
   real(DP), parameter   :: BETA = 2.85_DP !slope of sfd for remnants from LS12 2.85
   real(DP), parameter   :: C1 = 2.43_DP !LS12 constants
   real(DP), parameter   :: C2 = -0.0408_DP !LS12 constants
   real(DP), parameter   :: C3 = 1.86_DP !LS12 constants
   real(DP), parameter   :: C4 = 1.08_DP !LS12 constants
   real(DP), parameter   :: C5 = 2.5_DP !LS12 constants
   real(DP), parameter   :: CRUFU = 2.0_DP - 3 * MU_BAR ! central potential variable from Rufu and Aharonson (2019)
   real(DP), parameter   :: SUPERCAT_QRATIO = 1.8_DP ! See Section 4.1 of LS12
! Internals
   real(DP)           :: a1, alpha, aint, b, bcrit, egy, fgamma, l, lint, mu, phi, theta
   real(DP)           :: Qr, Qrd_pstar, Qr_erosion, Qr_supercat
   real(DP)           :: vcr, verosion, vescp, vhill, vimp, vsupercat
   real(DP)           :: mint, mtot
   real(DP)           :: rp, rhill 
   real(DP)           :: mresidual

! Executable code
   vimp = norm2(vb2(:) - vb1(:))
   b = calc_b(xh2, vb2, xh1, vb1)
   l = (rad1 + rad2) * (1 - b)
   egy = norm2(vb1)**2 / 2.0_DP - GC * Mcb / norm2(xh1)
   a1 = - GC * Mcb / 2.0_DP / egy
   mtot = m1 + m2 
   mu = (m1 * m2) / mtot
   if (l < 2 * rad2) then
      !calculate mint
      phi = 2 * acos((l - rad2) / rad2)
      aint = rad2**2 * (PI - (phi - sin(phi)) / 2.0_DP)
      lint = 2 * sqrt(rad2**2 - (rad2 - l / 2.0_DP) ** 2) 
      mint = aint * lint  ![kg]
      alpha = (l**2) * (3 * rad2 - l) / (4 * (rad2**3))
   else
      alpha = 1.0_DP
      mint = m2
   end if 
   rp = (3 * (m1 / den1 + alpha * m2 / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! (Mustill et al. 2018)
   !calculate vescp
   vescp = sqrt(2 * GC * mtot / rp) !Mustill et al. 2018 eq 6 
   !calculate rhill
   rhill = a1 * (m1 / 3.0_DP / (Mcb + m1))**(1.0_DP/3.0_DP)
   !calculate vhill
   if ((rad2 + rad1) < rhill) then 
     vhill = sqrt(2 * GC * m1 * ((rhill**2 - rhill * (rad1 + rad2)) / &
     (rhill**2 - 0.5_DP * (rad1 + rad2)**2)) / (rad1 + rad2))
   else
     vhill = vescp
   end if 
   !calculate Qr_pstar
   Qrd_pstar = calc_Qrd_pstar(m1, m2, alpha) * (vhill / vescp)**CRUFU !Rufu and Aharaonson eq (3)
   !calculate verosion
   Qr_erosion = 2 * (1.0_DP - m1 / mtot) * Qrd_pstar
   verosion = (2 * Qr_erosion * mtot / mu)** (1.0_DP / 2.0_DP)
   Qr = mu*(vimp**2) / mtot / 2.0_DP
   !calculate mass largest remnant mlr 
   mlr = (1.0_DP - Qr / Qrd_pstar / 2.0_DP) * mtot  ! [kg] # LS12 eq (5)
   !calculate vsupercat
   Qr_supercat = SUPERCAT_QRATIO * Qrd_pstar ! See LS12 Section 4.1 
   vsupercat = sqrt(2 * Qr_supercat * mtot / mu)
   !calculate vcr
   fgamma = (m1 - m2) / mtot
   theta = 1.0_DP - b
   vcr = vescp * (C1 * fgamma * theta**C5 + C2 * fgamma + C3 * theta**C5 + C4)
   bcrit = rad1 / (rad1 + rad2)

   if ((m1 < mtiny).or.(m2 < mtiny)) then 
     regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
     mlr = mtot
     mslr = 0.0_DP
     write(*,*) "FORCE MERGE"
   else 
      if( vimp < vescp) then
         regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
         mlr = mtot
         mslr = 0.0_DP
      else if (vimp < verosion) then 
         if (b < bcrit) then
            regime = COLLRESOLVE_REGIME_MERGE !partial accretion regime"
            mlr = mtot
            mslr = 0.0_DP
         else if ((b > bcrit) .and. (vimp < vcr)) then
            regime = COLLRESOLVE_REGIME_MERGE ! graze and merge
            mlr = mtot
            mslr = 0.0_DP
         else
            mlr = m1
            mslr = calc_Qrd_rev(m2,m1,mint,den1,den2,vimp)
            regime = COLLRESOLVE_REGIME_HIT_AND_RUN !hit and run
         end if 
      else if (vimp > verosion .and. vimp < vsupercat) then
         if (m2 < 0.001_DP * m1) then 
            regime = COLLRESOLVE_REGIME_MERGE !cratering regime"
            mlr = mtot
            mslr = 0.0_DP
         else 
            mslr = mtot * (3.0_DP - BETA) * (1.0_DP - N1 * mlr / mtot) / (N2 * BETA)  ! LS12 eq (37)
            regime = COLLRESOLVE_REGIME_DISRUPTION !disruption
         end if 
      else if (vimp > vsupercat) then 
         mlr = mtot * 0.1_DP * (Qr / (Qrd_pstar * SUPERCAT_QRATIO))**(-1.5_DP)   !LS12 eq (44)
         mslr = mtot * (3.0_DP - BETA) * (1.0_DP - N1 * mlr / mtot) / (N2 * BETA)  !LS12 eq (37)
         regime = COLLRESOLVE_REGIME_SUPERCATASTROPHIC ! supercatastrophic
      else 
         write(*,*) "Error no regime found in symba_regime"
      end if 
   end if 
   mresidual = mtot - mlr - mslr
   if (mresidual < 0.0_DP) then ! prevents final masses from going negative
      mlr = mlr + mresidual
   end if
      
   return 

! Internal functions
contains
   function calc_Qrd_pstar(mtarg,mp,alpha) result(Qrd_pstar)
      !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
      !!
      !! Calculates the corrected Q* for oblique impacts. See Eq. (15) of LS12.
      !!       Reference:
      !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
      !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
      !! 
      implicit none
      ! Arguments
      real(DP),intent(in) :: mtarg, mp, alpha
      ! Result
      real(DP)      :: Qrd_pstar
      ! Internals
      real(DP)      :: Qrd_star1, mu_alpha, mu, Qrd_star

      ! calc mu, mu_alpha
      mu = (mtarg * mp) / (mtarg + mp)  ! [kg]
      mu_alpha = (mtarg * alpha * mp) / (mtarg + alpha * mp)  ! [kg]
      ! calc Qrd_star1
      Qrd_star1 = (C_STAR * 4 * PI * DENSITY1 * GC * rp**2) / 5.0_DP
      ! calc Qrd_star
      Qrd_star = Qrd_star1 * (((mp / mtarg + 1.0_DP)**2) / (4 * mp / mtarg))**(2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP)  !(eq 23)
      ! calc Qrd_pstar, v_pstar
      Qrd_pstar = ((mu / mu_alpha)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP)) * Qrd_star  ! (eq 15)

      return
   end function calc_Qrd_pstar

   function calc_Qrd_rev(mp,mtarg,mint,den1,den2, vimp) result(mslr)
      !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
      !!
      !! Calculates mass of second largest fragment.
      !! 
      implicit none
      ! Arguments
      real(DP),intent(in) :: mp, mtarg, mint, den1, den2, vimp
      ! Result
      real(DP) :: mslr
      ! Internals
      real(DP) :: mtot_rev, mu_rev, gamma_rev, Qrd_star1, Qrd_star, mu_alpha_rev
      real(DP) :: Qrd_pstar, rC1, Qr_rev, Qrd_pstar_rev, Qr_supercat_rev

      ! calc mtlr, rC1, mu, gammalr
      mtot_rev =  mint + mp
      rC1 = (3 * (mint / den1 + mp / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! [m] Mustill et al 2018
      mu_rev = (mint * mp) / mtot_rev ! [kg] eq 49 LS12
      mu_alpha_rev = (mtarg * alpha * mp) / (mtarg + alpha * mp)
      gamma_rev = mint / mp ! eq 50 LS12
      !calc Qr_rev
      Qr_rev = mu_rev * (vimp**2) / (2 * mtot_rev)
      ! calc Qrd_star1, v_star1
      Qrd_star1 = (C_STAR * 4 * PI * mtot_rev * GC ) / rC1 / 5.0_DP
      ! calc Qrd_pstar_rev
      Qrd_star = Qrd_star1 * (((gamma_rev + 1.0_DP)**2) / (4 * gamma_rev)) ** (2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP) !(eq 52)
      Qrd_pstar = Qrd_star * ((mu_rev / mu_alpha_rev)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP))
      Qrd_pstar_rev = Qrd_pstar * (vhill / vescp)**CRUFU !Rufu and Aharaonson eq (3)
      !calc Qr_supercat_rev
      Qr_supercat_rev = 1.8_DP * Qrd_pstar_rev 
      if (Qr_rev > Qr_supercat_rev ) then 
         mslr = mtot_rev * (0.1_DP * ((Qr_rev / (Qrd_pstar_rev * 1.8_DP))**(-1.5_DP)))   !eq (44)
      else if ( Qr_rev < Qrd_pstar_rev ) then 
         mslr = mp 
      else 
         mslr = (1.0_DP - Qr_rev / Qrd_pstar_rev / 2.0_DP) * (mtot_rev)  ! [kg] #(eq 5)
      end if 

      if ( mslr > mp ) mslr = mp !check conservation of mass

      return
   end function calc_Qrd_rev

function calc_b(proj_pos, proj_vel, targ_pos, targ_vel) result(sintheta)
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Calculates the impact factor b = sin(theta), where theta is the angle between the relative velocity
   !! and distance vectors of the target and projectile bodies. See Fig. 2 of Leinhardt and Stewart (2012)
   !! 
   implicit none
   !! Arguments
   real(DP), dimension(:), intent(in) :: proj_pos, proj_vel, targ_pos, targ_vel
   !! Result
   real(DP)             :: sintheta
   !! Internals
   real(DP), dimension(NDIM)     :: imp_vel, distance, x_cross_v      

   imp_vel(:) = proj_vel(:) - targ_vel(:)
   distance(:) = proj_pos(:) - targ_pos(:)
   call util_crossproduct(distance, imp_vel, x_cross_v)
   sintheta = norm2(x_cross_v(:)) / norm2(distance(:)) / norm2(imp_vel(:))
   return 
end function calc_b

end subroutine symba_regime