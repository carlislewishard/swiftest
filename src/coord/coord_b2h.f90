subroutine coord_b2h(npl, swiftest_plA)
   !! Author: David Minton
   !! 
   !! Convert from barycentric to heliocentric coordinates, planets only
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine coord_b2h.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_b2h.f
! Modules
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => coord_b2h
   implicit none

! Arguments
   integer(I4B), intent(in)  :: npl
   type(swiftest_pl),intent(inout) :: swiftest_plA

! Internals
   integer(I4B)          :: i

! Executable code
   do i = 1, NDIM
      swiftest_plA%xh(i,1:npl) = swiftest_plA%xb(i,1:npl) - swiftest_plA%xb(i,1)
      swiftest_plA%vh(i,1:npl) = swiftest_plA%vb(i,1:npl) - swiftest_plA%vb(i,1)
   end do

   return

end subroutine coord_b2h