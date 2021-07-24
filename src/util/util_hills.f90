subroutine util_hills(npl, swiftest_plA)
   !! Author: David A. Minton
   !!
   !! Compute Hill sphere radii of planet. in the case of hyperbolic orbits, use the heliocentric radius instead of semimajor axis
   !!
   !! Adapted from David E. Kaufmann's Swifter routine util_hills.f90
   !! Adapted from Hal Levison's Swift routine util_hills.f
! Modules
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => util_hills
   implicit none

! arguments
   integer(I4B), intent(in)  :: npl
   type(swiftest_pl), intent(inout) :: swiftest_plA

! internals
   integer(I4B)        :: i
   real(DP)            :: mu, energy, ap, r, v2, hx, hy, hz, h2, fac, ecc

! executable code
   associate(GMcb => swiftest_plA%mass(1), GMp => swiftest_plA%mass, xh => swiftest_plA%xh, vh => swiftest_plA%vh, rhill => swiftest_plA%rhill)
      do i = 2, npl
         if (GMp(i) > 0.0_DP) then
            if ((swiftest_plA%status(i) == DISRUPTION) .or. (swiftest_plA%status(i) == SUPERCATASTROPHIC) .or. (swiftest_plA%status(i) == HIT_AND_RUN)) then
               hx = xh(2,i) * vh(3,i) - xh(3,i) * vh(2,i)
               hy = xh(3,i) * vh(1,i) - xh(1,i) * vh(3,i)   
               hz = xh(1,i) * vh(2,i) - xh(2,i) * vh(1,i)
               h2 = hx*hx + hy*hy + hz*hz
               mu = GMcb + GMp(i)
               r = norm2(xh(:, i))
               v2 = dot_product(vh(:, i), vh(:, i))
               energy = 0.5_DP * v2 - mu / r
               if (energy < 0.0_DP) then
                  ap = -0.5_DP * mu / energy
                  fac = 1.0_DP - h2/(mu*ap)
                  ecc = sqrt(fac)
                  if (ecc > 0.5_DP) then
                     ap = r ! use the heliocentric radius for highly eccentric bodies
                  end if
               else
                  ap = r ! use the heliocentric radius for the hill radius of hyperbolic orbits 
                       !(this is probably good enough for most purposes, but probably worth investigating)
               end if
               rhill(i) = ap * (((GMp(i) / mu) / 3.0_DP)**(1.0_DP / 3.0_DP))
            end if
         else
            rhill(i) = 0.0_DP
         end if
      end do
   end associate

   return

end subroutine util_hills