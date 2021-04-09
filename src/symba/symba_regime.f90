subroutine symba_regime(Mcb, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, mtiny, Qresidual)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Determine the collisional regime of two colliding bodies. 
   !! Current version requires all values to be converted to SI units prior to calling the function
   !!       References:
   !!       Kokubo, E., Genda, H., 2010. Formation of Terrestrial Planets from Protoplanets Under a Realistic Accretion 
   !!          Condition. ApJL 714, L21. https://doi.org/10.1088/2041-8205/714/1/L21
   !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
   !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
   !!       Mustill, A.J., Davies, M.B., Johansen, A., 2018. The dynamical evolution of transiting planetary systems including 
   !!          a realistic collision prescription. Mon Not R Astron Soc 478, 2896–2908. https://doi.org/10.1093/mnras/sty1273
   !!       Rufu, R., Aharonson, O., 2019. Impact Dynamics of Moons Within a Planetary Potential. J. Geophys. Res. Planets 124, 
   !!          1008–1019. https://doi.org/10.1029/2018JE005798
   !!       Stewart, S.T., Leinhardt, Z.M., 2012. Collisions between Gravity-dominated Bodies. II. The Diversity of Impact 
   !!          Outcomes during the End Stage of Planet Formation. ApJ 751, 32. https://doi.org/10.1088/0004-637X/751/1/32
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
   real(DP), intent(out)          :: Mlr, Mslr
   real(DP), intent(in)           :: Mcb, m1, m2, rad1, rad2, den1, den2, mtiny 
   real(DP), dimension(:), intent(in)   :: xh1, xh2, vb1, vb2
   real(DP), intent(out)          :: Qresidual !! The residual energy after the collision 
! Constants
   integer(I4B), parameter :: N1 = 1  !number of objects with mass equal to the largest remnant from LS12
   integer(I4B), parameter :: N2 = 2  !number of objects with mass larger than second largest remnant from LS12
   real(DP), parameter   :: DENSITY1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
   real(DP), parameter   :: MU_BAR = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
   real(DP), parameter   :: BETA = 2.85_DP !slope of sfd for remnants from LS12 2.85
   real(DP), parameter   :: C1 = 2.43_DP  !! Kokubo & Genda (2010) eq. (3)
   real(DP), parameter   :: C2 = -0.0408_DP !! Kokubo & Genda (2010) eq. (3)
   real(DP), parameter   :: C3 = 1.86_DP !! Kokubo & Genda (2010) eq. (3)
   real(DP), parameter   :: C4 = 1.08_DP !! Kokubo & Genda (2010) eq. (3)
   real(DP), parameter   :: CRUFU = 2.0_DP - 3 * MU_BAR ! central potential variable from Rufu and Aharonson (2019)
   real(DP), parameter   :: SUPERCAT_QRATIO = 1.8_DP ! See Section 4.1 of LS12
! Internals
   real(DP)           :: a1, alpha, aint, b, bcrit, c_star, egy, zeta, l, lint, mu, phi, theta
   real(DP)           :: Qr, Qrd_pstar, Qr_erosion, Qr_supercat
   real(DP)           :: Vhr, Verosion, Vescp, Vhill, Vimp, Vsupercat
   real(DP)           :: Mint, Mtot
   real(DP)           :: Rp, rhill 
   real(DP)           :: Mresidual

! Executable code
   Vimp = norm2(vb2(:) - vb1(:))
   b = calc_b(xh2, vb2, xh1, vb1)
   l = (rad1 + rad2) * (1 - b)
   egy = 0.5_DP * dot_product(vb1, vb1) - GC * Mcb / norm2(xh1)
   a1 = - GC * Mcb / 2.0_DP / egy
   Mtot = m1 + m2 
   mu = (m1 * m2) / Mtot
   if (l < 2 * rad2) then
      !calculate Mint
      phi = 2 * acos((l - rad2) / rad2)
      aint = rad2**2 * (PI - (phi - sin(phi)) / 2.0_DP)
      lint = 2 * sqrt(rad2**2 - (rad2 - l / 2.0_DP) ** 2) 
      Mint = aint * lint  ![kg]
      alpha = (l**2) * (3 * rad2 - l) / (4 * (rad2**3))
   else
      alpha = 1.0_DP
      Mint = m2
   end if 
   Rp = (3 * (m1 / den1 + alpha * m2 / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! (Mustill et al. 2018)
   c_star = calc_c_star(Rp)
   !calculate Vescp
   Vescp = sqrt(2 * GC * Mtot / Rp) !Mustill et al. 2018 eq 6 
   !calculate rhill
   rhill = a1 * (m1 / 3.0_DP / (Mcb + m1))**(1.0_DP/3.0_DP)
   !calculate Vhill
   if ((rad2 + rad1) < rhill) then 
     Vhill = sqrt(2 * GC * m1 * ((rhill**2 - rhill * (rad1 + rad2)) / &
     (rhill**2 - 0.5_DP * (rad1 + rad2)**2)) / (rad1 + rad2))
   else
     Vhill = Vescp
   end if 
   !calculate Qr_pstar
   Qrd_pstar = calc_Qrd_pstar(m1, m2, alpha, c_star) * (Vhill / Vescp)**CRUFU !Rufu and Aharaonson eq (3)
   !calculate Verosion
   Qr_erosion = 2 * (1.0_DP - m1 / Mtot) * Qrd_pstar
   Verosion = (2 * Qr_erosion * Mtot / mu)** (1.0_DP / 2.0_DP)
   Qr = mu*(Vimp**2) / Mtot / 2.0_DP
   !calculate mass largest remnant Mlr 
   Mlr = (1.0_DP - Qr / Qrd_pstar / 2.0_DP) * Mtot  ! [kg] # LS12 eq (5)
   !calculate Vsupercat
   Qr_supercat = SUPERCAT_QRATIO * Qrd_pstar ! See LS12 Section 4.1 
   Vsupercat = sqrt(2 * Qr_supercat * Mtot / mu)
   !calculate Vhr
   zeta = (m1 - m2) / Mtot
   theta = 1.0_DP - b
   Vhr = Vescp * (C1 * zeta**2 * theta**(2.5_DP) + C2 * zeta**2 + C3 * theta**(2.5_DP) + C4) ! Kokubo & Genda (2010) eq. (3)
   bcrit = rad1 / (rad1 + rad2)

   if ((m1 < mtiny).or.(m2 < mtiny)) then 
     regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
     Mlr = Mtot
     Mslr = 0.0_DP
     write(*,*) "FORCE MERGE"
   else 
      if( Vimp < Vescp) then
         regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
         Mlr = Mtot
         Mslr = 0.0_DP
         Qresidual = 0.0_DP
      else if (Vimp < Verosion) then 
         if (b < bcrit) then
            regime = COLLRESOLVE_REGIME_MERGE !partial accretion regime"
            Mlr = Mtot
            Mslr = 0.0_DP
            Qresidual = 0.0_DP
         else if ((b > bcrit) .and. (Vimp < Vhr)) then
            regime = COLLRESOLVE_REGIME_MERGE ! graze and merge
            Mlr = Mtot
            Mslr = 0.0_DP
            Qresidual = 0.0_DP
         else
            Mlr = m1
            Mslr = calc_Qrd_rev(m2, m1, Mint, den1, den2, Vimp, c_star)
            regime = COLLRESOLVE_REGIME_HIT_AND_RUN !hit and run
            Qresidual = Qr
         end if 
      else if (Vimp > Verosion .and. Vimp < Vsupercat) then
         if (m2 < 0.001_DP * m1) then 
            regime = COLLRESOLVE_REGIME_MERGE !cratering regime"
            Mlr = Mtot
            Mslr = 0.0_DP
            Qresidual = 0.0_DP
         else 
            Mslr = Mtot * (3.0_DP - BETA) * (1.0_DP - N1 * Mlr / Mtot) / (N2 * BETA)  ! LS12 eq (37)
            regime = COLLRESOLVE_REGIME_DISRUPTION !disruption
            Qresidual = Qr - Qr_erosion
         end if 
      else if (Vimp > Vsupercat) then 
         Mlr = Mtot * 0.1_DP * (Qr / (Qrd_pstar * SUPERCAT_QRATIO))**(-1.5_DP)   !LS12 eq (44)
         Mslr = Mtot * (3.0_DP - BETA) * (1.0_DP - N1 * Mlr / Mtot) / (N2 * BETA)  !LS12 eq (37)
         regime = COLLRESOLVE_REGIME_SUPERCATASTROPHIC ! supercatastrophic
         Qresidual = Qr - Qr_supercat
      else 
         write(*,*) "Error no regime found in symba_regime"
      end if 
   end if 
   Mresidual = Mtot - Mlr - Mslr
   if (Mresidual < 0.0_DP) then ! prevents final masses from going negative
      Mlr = Mlr + Mresidual
   end if
      
   return 

! Internal functions
contains
   function calc_Qrd_pstar(Mtarg, Mp, alpha, c_star) result(Qrd_pstar)
      !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
      !!
      !! Calculates the corrected Q* for oblique impacts. See Eq. (15) of LS12.
      !!       Reference:
      !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
      !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
      !! 
      implicit none
      ! Arguments
      real(DP),intent(in) :: Mtarg, Mp, alpha, c_star
      ! Result
      real(DP)      :: Qrd_pstar
      ! Internals
      real(DP)      :: Qrd_star1, mu_alpha, mu, Qrd_star

      ! calc mu, mu_alpha
      mu = (Mtarg * Mp) / (Mtarg + Mp)  ! [kg]
      mu_alpha = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)  ! [kg]
      ! calc Qrd_star1
      Qrd_star1 = (c_star * 4 * PI * DENSITY1 * GC * Rp**2) / 5.0_DP
      ! calc Qrd_star
      Qrd_star = Qrd_star1 * (((Mp / Mtarg + 1.0_DP)**2) / (4 * Mp / Mtarg))**(2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP)  !(eq 23)
      ! calc Qrd_pstar, v_pstar
      Qrd_pstar = ((mu / mu_alpha)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP)) * Qrd_star  ! (eq 15)

      return
   end function calc_Qrd_pstar

   function calc_Qrd_rev(Mp, Mtarg, Mint, den1, den2, Vimp, c_star) result(Mslr)
      !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
      !!
      !! Calculates mass of second largest fragment.
      !! 
      implicit none
      ! Arguments
      real(DP),intent(in) :: Mp, Mtarg, Mint, den1, den2, Vimp, c_star
      ! Result
      real(DP) :: Mslr
      ! Internals
      real(DP) :: mtot_rev, mu_rev, gamma_rev, Qrd_star1, Qrd_star, mu_alpha_rev
      real(DP) :: Qrd_pstar, Rc1, Qr_rev, Qrd_pstar_rev, Qr_supercat_rev

      ! calc Mslr, Rc1, mu, gammalr
      mtot_rev =  Mint + Mp
      Rc1 = (3 * (Mint / den1 + Mp / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! [m] Mustill et al 2018
      mu_rev = (Mint * Mp) / mtot_rev ! [kg] eq 49 LS12
      mu_alpha_rev = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)
      gamma_rev = Mint / Mp ! eq 50 LS12
      !calc Qr_rev
      Qr_rev = mu_rev * (Vimp**2) / (2 * mtot_rev)
      ! calc Qrd_star1, v_star1
      Qrd_star1 = (c_star * 4 * PI * mtot_rev * GC ) / Rc1 / 5.0_DP
      ! calc Qrd_pstar_rev
      Qrd_star = Qrd_star1 * (((gamma_rev + 1.0_DP)**2) / (4 * gamma_rev)) ** (2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP) !(eq 52)
      Qrd_pstar = Qrd_star * ((mu_rev / mu_alpha_rev)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP))
      Qrd_pstar_rev = Qrd_pstar * (Vhill / Vescp)**CRUFU !Rufu and Aharaonson eq (3)
      !calc Qr_supercat_rev
      Qr_supercat_rev = 1.8_DP * Qrd_pstar_rev 
      if (Qr_rev > Qr_supercat_rev ) then 
         Mslr = mtot_rev * (0.1_DP * ((Qr_rev / (Qrd_pstar_rev * 1.8_DP))**(-1.5_DP)))   !eq (44)
      else if ( Qr_rev < Qrd_pstar_rev ) then 
         Mslr = Mp 
      else 
         Mslr = (1.0_DP - Qr_rev / Qrd_pstar_rev / 2.0_DP) * (mtot_rev)  ! [kg] #(eq 5)
      end if 

      if ( Mslr > Mp ) Mslr = Mp !check conservation of mass

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

   function calc_c_star(Rc1) result(c_star)
      !! author: David A. Minton
      !!
      !! Calculates c_star as a function of impact equivalent radius. It inteRpolates between 5 for ~1 km sized bodies to
      !! 1.8 for ~10000 km sized bodies. See LS12 Fig. 4 for details.
      !! 
      implicit none
      !! Arguments
      real(DP), intent(in) :: Rc1
      !! Result
      real(DP)             :: c_star
      !! Internals
      real(DP), parameter  :: loR   = 1.0e3_DP ! Lower bound of inteRpolation size (m)
      real(DP), parameter  :: hiR   = 1.0e7_DP ! Upper bound of inteRpolation size (m)
      real(DP), parameter  :: loval = 5.0_DP   ! Value of C* at lower bound
      real(DP), parameter  :: hival = 1.8_DP   ! Value of C* at upper bound

      if (Rc1 < loR) then
         c_star = loval
      else if (Rc1 < hiR) then
         c_star = loval + (hival - loval) * log(Rc1 / loR) / log(hiR /loR)
      else
         c_star = hival
      end if
      return
   end function calc_c_star
end subroutine symba_regime