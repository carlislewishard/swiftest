submodule (swiftest_classes) s_discard
   use swiftest
contains

   module subroutine discard_system(self, param)
      !! author: David A. Minton
      !!
      !! Calls the discard methods for each body class and then the write method if any discards were detected
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical :: lany_discards

      associate(system => self, tp => self%tp, pl => self%pl)
         lany_discards = .false.
         call pl%discard(system, param)
         call tp%discard(system, param)
         if (tp%nbody > 0) lany_discards = lany_discards .or.  any(tp%ldiscard(:))
         if (pl%nbody > 0) lany_discards = lany_discards .or.  any(pl%ldiscard(:))
         if (lany_discards) call system%write_discard(param)
      end associate

      return
   end subroutine discard_system


   module subroutine discard_pl(self, system, param)
      !! author: David A. Minton
      !!
      !!  Placeholder method for discarding massive bodies. This method does nothing except to ensure that the discard flag is set to false. 
      !!  This method is intended to be overridden by more advanced integrators.
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameter

      self%ldiscard(:) = .false.

      return
   end subroutine discard_pl


   module subroutine discard_tp(self, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if particles should be discarded based on their positions relative to the massive bodies
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
      !! Adapted from Hal Levison's Swift routine discard.
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameter

      associate(tp => self, ntp => self%nbody, cb => system%cb, pl => system%pl, npl => system%pl%nbody)
         if (ntp == 0) return 
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
             (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
            if (npl > 0) call pl%h2b(cb) 
            if (ntp > 0) call tp%h2b(cb) 
         end if
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) then
            if (ntp > 0) call discard_sun_tp(tp, system, param)
         end if
         if (param%qmin >= 0.0_DP .and. ntp > 0) call discard_peri_tp(tp, system, param)
         if (param%lclose .and. ntp > 0) call discard_pl_tp(tp, system, param)
         if (any(tp%ldiscard)) call tp%spill(system%tp_discards, tp%ldiscard)
      end associate

      return
   end subroutine discard_tp


   subroutine discard_sun_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !!  Check to see if test particles should be discarded based on their positions relative to the Sun
      !!        or because they are unbound from the system
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_sun.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)        :: i
      real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2

      associate(ntp => tp%nbody, cb => system%cb, t => param%t, msys => system%msys)
         rmin2 = max(param%rmin * param%rmin, cb%radius * cb%radius)
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               rh2 = dot_product(tp%xh(:, i), tp%xh(:, i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  tp%status(i) = DISCARDED_RMAX
                  write(*, *) "Particle ", tp%id(i), " too far from sun at t = ", t
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  tp%status(i) = DISCARDED_RMIN
                  write(*, *) "Particle ", tp%id(i), " too close to sun at t = ", t
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(tp%xb(:, i),  tp%xb(:, i))
                  vb2 = dot_product(tp%vb(:, i), tp%vb(:, i))
                  energy = 0.5_DP * vb2 - msys / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     tp%status(i) = DISCARDED_RMAXU
                     write(*, *) "Particle ", tp%id(i), " is unbound and too far from barycenter at t = ", t
                     tp%ldiscard(i) = .true.
                     tp%lmask(i) = .false.
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine discard_sun_tp


   subroutine discard_peri_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameterss
      ! Internals
      logical, save             :: lfirst = .true.
      integer(I4B)              :: i, j, ih
      real(DP)                  :: r2
      real(DP), dimension(NDIM) :: dx
   
      associate(cb => system%cb, ntp => tp%nbody, pl => system%pl, npl => system%pl%nbody, t => param%t)
         call tp%get_peri(system, param)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               if (tp%isperi(i) == 0) then
                  ih = 1
                  do j = 1, npl
                     dx(:) = tp%xh(:, i) - pl%xh(:, j)
                     r2 = dot_product(dx(:), dx(:))
                     if (r2 <= (pl%rhill(j))**2) ih = 0
                  end do
                  if (ih == 1) then
                     if ((tp%atp(i) >= param%qmin_alo) .and.    &
                        (tp%atp(i) <= param%qmin_ahi) .and.    &           
                        (tp%peri(i) <= param%qmin)) then
                        tp%status(i) = DISCARDED_PERI
                        write(*, *) "Particle ", tp%id(i), " perihelion distance too small at t = ", t
                        tp%ldiscard(i) = .true.
                     end if
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine discard_peri_tp


   subroutine discard_pl_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals 
      integer(I4B)              :: i, j, isp
      real(DP)                  :: r2min, radius
      real(DP), dimension(NDIM) :: dx, dv
   
      associate(ntp => tp%nbody, pl => system%pl, npl => system%pl%nbody, t => param%t, dt => param%dt)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%xh(:, i) - pl%xh(:, j)
                  dv(:) = tp%vh(:, i) - pl%vh(:, j)
                  radius = pl%radius(j)
                  call discard_pl_close(dx(:), dv(:), dt, radius**2, isp, r2min)
                  if (isp /= 0) then
                     tp%status(i) = DISCARDED_PLR
                     tp%lmask(i) = .false.
                     pl%ldiscard(j) = .true.
                     write(*, *) "Particle ", tp%id(i), " too close to massive body ", pl%id(j), " at t = ", t
                     tp%ldiscard(i) = .true.
                     exit
                  end if
               end do
            end if
         end do
      end associate

      return
   end subroutine discard_pl_tp
   

   subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
      !! author: David A. Minton
      !!
      !!  Check to see if a test particle and massive body are having, or will have within the next time step, an encounter such
      !!          that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl_close.f90
      !! Adapted from Hal Levison's Swift routine discard_pl_close.f
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)    :: dx, dv
      real(DP), intent(in)                  :: dt, r2crit
      integer(I4B), intent(out)             :: iflag
      real(DP), intent(out)                 :: r2min
      ! Internals
      real(DP) :: r2, v2, vdotr, tmin
      
      r2 = dot_product(dx(:), dx(:))
      if (r2 <= r2crit) then
         iflag = 1
      else
         vdotr = dot_product(dx(:), dv(:))
         if (vdotr > 0.0_DP) then
            iflag = 0
         else
            v2 = dot_product(dv(:), dv(:))
            tmin = -vdotr / v2
            if (tmin < dt) then
               r2min = r2 - vdotr * vdotr / v2
            else
               r2min = r2 + 2 * vdotr * dt + v2 * dt**2
            end if
            r2min = min(r2min, r2)
            if (r2min <= r2crit) then
               iflag = 1
            else
               iflag = 0
            end if
         end if
      end if
   
      return
   end subroutine discard_pl_close

end submodule s_discard
