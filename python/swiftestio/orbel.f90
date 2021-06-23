! Created by  on 6/22/21.

module orbel
   implicit none

contains

   subroutine xv2el(mu, x, v, a, e, inc, capom, omega, capm)
      !! author: David A. Minton
      !!
      !! Compute osculating orbital elements from relative Cartesian position and velocity
      !!  All angular measures are returned in radians
      !!      If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
      !!
      !!      If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
      !!
      !!      References: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 201 - 206.
      !!              Fitzpatrick, P. M. 1970. Principles of Celestial Mechanics, (Academic Press), 69 - 73.
      !!              Roy, A. E. 1982. Orbital Motion, (Adam Hilger, Ltd.), 75 - 95
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2el.f90
      !! Adapted from Martin Duncan's Swift routine orbel_xv2el.f
      implicit none
      ! Arguments
      real*8, intent(in)  :: mu
      real*8, dimension(:), intent(in)  :: x, v
      real*8, intent(out) :: a, e, inc, capom, omega, capm
      !f2py intent(in) mu
      !f2py intent(in) x
      !f2py intent(in) v
      !f2py intent(out) a
      !f2py intent(out) e
      !f2py intent(out) inc
      !f2py intent(out) capom
      !f2py intent(out) omega
      !f2py intent(out) capm

      ! Internals

      integer :: iorbit_type
      real*8   :: r, v2, h2, h, rdotv, energy, fac, u, w, cw, sw, face, cape, tmpf, capf
      real*8, dimension(3) :: hvec
      integer, parameter :: ELLIPSE   = -1 !! Symbolic names for orbit types - ellipse
      integer, parameter :: PARABOLA  =  0 !! Symbolic names for orbit types - parabola
      integer, parameter :: HYPERBOLA =  1 !! Symbolic names for orbit types - hyperbola
      integer, parameter :: DP = kind(1.d0)
      real*8,     parameter :: VSMALL    = 4.0E-15_DP
      real*8, parameter :: PI     = 3.141592653589793238462643383279502884197_DP !! Definition of /(\pi\)
      real*8, parameter :: TWOPI  = 6.283185307179586476925286766559005768394_DP !! Definition of 2 \pi

      a = 0.0_DP
      e = 0.0_DP
      inc = 0.0_DP
      capom = 0.0_DP
      omega = 0.0_DP
      capm = 0.0_DP
      r = sqrt(dot_product(x(:), x(:)))
      v2 = dot_product(v(:), v(:))
      hvec = crossproduct(x(:), v(:))
      h2 = dot_product(hvec(:), hvec(:))
      h = sqrt(h2)
      if (h2 == 0.0_DP) return
      rdotv = dot_product(x(:), v(:))
      energy = 0.5_DP * v2 - mu / r
      fac = hvec(3) / h
      if (fac < -1.0_DP) then
         inc = PI
      else if (fac < 1.0_DP) then
         inc = acos(fac)
      end if
      fac = sqrt(hvec(1)**2 + hvec(2)**2) / h
      if (fac**2 < VSMALL) then
         u = atan2(x(2), x(1))
         if (hvec(3) < 0.0_DP) u = -u
      else
         capom = atan2(hvec(1), -hvec(2))
         u = atan2(x(3) / sin(inc), x(1) * cos(capom) + x(2) * sin(capom))
      end if
      if (capom < 0.0_DP) capom = capom + TWOPI
      if (u < 0.0_DP) u = u + TWOPI
      if (abs(energy * r / mu) < sqrt(VSMALL)) then
         iorbit_type = PARABOLA
      else
         a = -0.5_DP * mu / energy
         if (a < 0.0_DP) then
            fac = -h2 / (mu * a)
            if (fac > VSMALL) then
               iorbit_type = HYPERBOLA
            else
               iorbit_type = PARABOLA
            end if
         else
            iorbit_type = ELLIPSE
         end if
      end if
      select case (iorbit_type)
         case (ELLIPSE)
            fac = 1.0_DP - h2 / (mu * a)
            if (fac > VSMALL) then
               e = sqrt(fac)
               cape = 0.0_DP
               face = (a - r) / (a * e)
               if (face < -1.0_DP) then
                  cape = PI
               else if (face < 1.0_DP) then
                  cape = acos(face)
               end if
               if (rdotv < 0.0_DP) cape = TWOPI - cape
               fac = 1.0_DP - e * cos(cape)
               cw = (cos(cape) - e) / fac
               sw = sqrt(1.0_DP - e**2) * sin(cape) / fac
               w = atan2(sw, cw)
               if (w < 0.0_DP) w = w + TWOPI
            else
               cape = u
               w = u
            end if
            capm = cape - e * sin(cape)
         case (PARABOLA)
            a = 0.5_DP * h2 / mu
            e = 1.0_DP
            w = 0.0_DP
            fac = 2 * a / r - 1.0_DP
            if (fac < -1.0_DP) then
               w = PI
            else if (fac < 1.0_DP) then
               w = acos(fac)
            end if
            if (rdotv < 0.0_DP) w = TWOPI - w
            tmpf = tan(0.5_DP * w)
            capm = tmpf * (1.0_DP + tmpf * tmpf / 3.0_DP)
         case (HYPERBOLA)
            e = sqrt(1.0_DP + fac)
            tmpf = max((a - r) / (a * e), 1.0_DP)
            capf = log(tmpf + sqrt(tmpf**2 - 1.0_DP))
            if (rdotv < 0.0_DP) capf = -capf
            fac = e * cosh(capf) - 1.0_DP
            cw = (e - cosh(capf)) / fac
            sw = sqrt(e * e - 1.0_DP) * sinh(capf) / fac
            w = atan2(sw, cw)
            if (w < 0.0_DP) w = w + TWOPI
            capm = e * sinh(capf) - capf
      end select
      omega = u - w
      if (omega < 0.0_DP) omega = omega + TWOPI

      return
      contains
         function crossproduct(A, B) result(C)
            implicit none
            real*8, dimension(:), intent(in) :: A, B
            real*8, dimension(3) :: C
            C(1) = A(2) * B(3) - A(3) * B(2)
            C(2) = A(3) * B(1) - A(1) * B(3)
            C(3) = A(1) * B(2) - A(2) * B(1)
            return
         end function crossproduct
   end subroutine xv2el

end module orbel