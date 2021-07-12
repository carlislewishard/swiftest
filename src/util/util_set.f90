submodule(swiftest_classes) s_util_set
   !! author: David A. Minton
   !! This submodule contains a collection of setter method implementations
   use swiftest
contains

   module subroutine util_set_beg_end_cb(self, aoblbeg, aoblend)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of aoblbeg and aoblend 
      implicit none
      ! Arguments
      class(swiftest_cb),     intent(inout)          :: self    !! Swiftest central body object
      real(DP), dimension(:), intent(in),   optional :: aoblbeg !! Oblateness acceleration term at beginning of step
      real(DP), dimension(:), intent(in),   optional :: aoblend !! Oblateness acceleration term at end of step

      if (present(aoblbeg)) self%aoblbeg = aoblbeg
      if (present(aoblend)) self%aoblend = aoblend
      return

   end subroutine util_set_beg_end_cb
   
   module subroutine util_set_beg_end_pl(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
      implicit none
      ! Arguments
      class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
      real(DP), dimension(:,:), intent(in),   optional :: xbeg, xend, vbeg

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end subroutine util_set_beg_end_pl

   module subroutine util_set_msys(self)
      !! author: David A. Minton
      !!
      !! Sets the value of msys and the vector mass quantities based on the total mass of the system
      implicit none
      class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system objec
      self%msys = self%cb%mass + sum(self%pl%mass(1:self%pl%nbody))

      return
   end subroutine util_set_msys

   module subroutine util_set_mu_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      implicit none
      class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody > 0) self%mu(:) = cb%Gmass + self%Gmass(:)

      return
   end subroutine util_set_mu_pl

   module subroutine util_set_mu_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      implicit none
      class(swiftest_tp),           intent(inout) :: self !! Swiftest test particle object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody > 0) self%mu(:) = cb%Gmass

      return
   end subroutine util_set_mu_tp

   module subroutine util_set_rhill(self,cb)
      !! author: David A. Minton
      !!
      !! Sets the value of the Hill's radius
      implicit none
      class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest massive body object

      if (self%nbody > 0) then
         call self%xv2el(cb) 
         self%rhill(:) = self%a(:) * (self%Gmass(:) / cb%Gmass / 3)**THIRD 
      end if

      return
   end subroutine util_set_rhill

   module subroutine util_set_ir3h(self)
      !! author: David A. Minton
      !!
      !! Sets the inverse heliocentric radius term (1/rh**3) for all bodies in a structure
      implicit none
      class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object

      integer(I4B) :: i
      real(DP) :: r2, irh

      if (self%nbody > 0) then

         do i = 1, self%nbody
            r2 = dot_product(self%xh(:, i), self%xh(:, i))
            irh = 1.0_DP / sqrt(r2)
            self%ir3h(i) = irh / r2
         end do
      end if

      return
   end subroutine util_set_ir3h
end submodule s_util_set