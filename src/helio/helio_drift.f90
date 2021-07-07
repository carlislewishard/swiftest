submodule (helio_classes) s_helio_drift
   use swiftest
contains
   module subroutine helio_drift_pl(self, system, param, dt)

      !! author: David A. Minton
      !!
      !! Loop through massive bodies and call Danby drift routine
      !! New vectorized version included
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_drift.f90
      !! Adapted from Hal Levison's Swift routine drift.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: dt     !! Stepsize)
      ! Internals
      integer(I4B) :: i !! Loop counter
      real(DP) :: rmag, vmag2, energy
      integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
      real(DP), dimension(:), allocatable          :: dtp

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return

         allocate(iflag(npl))
         iflag(:) = 0
         allocate(dtp(npl))

         if (param%lgr) then
            do i = 1,npl
               rmag = norm2(pl%xh(:, i))
               vmag2 = dot_product(pl%vb(:, i),  pl%vb(:, i))
               energy = 0.5_DP * vmag2 - pl%mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 

         call drift_one(pl%mu(1:npl), pl%xh(1,1:npl), pl%xh(2,1:npl), pl%xh(3,1:npl), &
                                      pl%vb(1,1:npl), pl%vb(2,1:npl), pl%vb(3,1:npl), &
                                      dtp(1:npl), iflag(1:npl))
         if (any(iflag(1:npl) /= 0)) then
            do i = 1, npl
               write(*, *) " Planet ", pl%name(i), " is lost!!!!!!!!!!"
               write(*, *) pl%xh(:,i)
               write(*, *) pl%vb(:,i)
               write(*, *) " stopping "
               call util_exit(FAILURE)
            end do
         end if
      end associate
      return
   end subroutine helio_drift_pl
   
   module subroutine helio_drift_linear_pl(self, system, dt, pt)
      !! author: David A. Minton
      !!
      !! Perform linear drift of massive bodies due to barycentric momentum of Sun
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift.f90
      !! Adapted from Hal Levison's Swift routine helio_lindrift.f
      implicit none
      ! Arguments
      class(helio_pl),           intent(inout) :: self   !! Helio massive body object
      class(helio_nbody_system), intent(in)    :: system !! Swiftest nbody system object
      real(DP),                  intent(in)    :: dt     !! Stepsize
      real(DP), dimension(:),    intent(out)   :: pt     !! negative barycentric velocity of the central body
  
      associate(pl => self, npl => self%nbody, cb => system%cb)
         pt(1) = sum(pl%Gmass(1:npl) * pl%vb(1,1:npl))
         pt(2) = sum(pl%Gmass(1:npl) * pl%vb(2,1:npl))
         pt(3) = sum(pl%Gmass(1:npl) * pl%vb(3,1:npl))
         pt(:) = pt(:) / cb%Gmass
         pl%xh(1,1:npl) = pl%xh(1,1:npl) + pt(1) * dt
         pl%xh(2,1:npl) = pl%xh(2,1:npl) + pt(2) * dt
         pl%xh(3,1:npl) = pl%xh(3,1:npl) + pt(3) * dt
      end associate
   
      return
   end subroutine helio_drift_linear_pl

   module subroutine helio_drift_tp(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop through test particles and call Danby drift routine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_drift_tp.f90
      !! Adapted from Hal Levison's Swift routine drift_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Helio test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B) :: i !! Loop counter
      real(DP) :: rmag, vmag2, energy
      real(DP), dimension(:), allocatable    :: dtp
      integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
   
      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         allocate(iflag(ntp))
         allocate(dtp(ntp))
         iflag(:) = 0

         if (param%lgr) then
            do i = 1,ntp
               rmag = norm2(tp%xh(:, i))
               vmag2 = dot_product(tp%vh(:, i), tp%vh(:, i))
               energy = 0.5_DP * vmag2 - tp%mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 
         call drift_one(tp%mu(1:ntp), tp%xh(1,1:ntp), tp%xh(2,1:ntp), tp%xh(3,1:ntp), &
                                      tp%vh(1,1:ntp), tp%vh(2,1:ntp), tp%vh(3,1:ntp), &
                                      dtp(1:ntp), iflag(1:ntp))
         if (any(iflag(1:ntp) /= 0)) then
            tp%status = DISCARDED_DRIFTERR
            do i = 1, ntp
               if (iflag(i) /= 0) write(*, *) "Particle ", tp%name(i), " lost due to error in Danby drift"
            end do
         end if
      end associate

      return
   end subroutine helio_drift_tp

   module subroutine helio_drift_linear_tp(self, system, dt, pt)
      !! author: David A. Minton
      !!
      !! Perform linear drift of test particles due to barycentric momentum of Sun
      !! New vectorized version included
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
      implicit none
      ! Arguments
      class(helio_tp),           intent(inout) :: self   !! Helio test particleb object
      class(helio_nbody_system), intent(in)    :: system !! Swiftest nbody system object
      real(DP),                  intent(in)    :: dt     !! Stepsize
      real(DP), dimension(:),    intent(in)    :: pt     !! negative barycentric velocity of the central body
  
      associate(tp => self, ntp => self%nbody)
         where (tp%status(1:ntp) == ACTIVE)
            tp%xh(1, 1:ntp) = tp%xh(1, 1:ntp) + pt(1) * dt
            tp%xh(2, 1:ntp) = tp%xh(2, 1:ntp) + pt(2) * dt
            tp%xh(3, 1:ntp) = tp%xh(3, 1:ntp) + pt(3) * dt
         end where
      end associate
   
      return
   end subroutine helio_drift_linear_tp

end submodule s_helio_drift
