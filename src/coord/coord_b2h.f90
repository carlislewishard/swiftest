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
   real(DP), dimension(NDIM) :: xtmp, vtmp

! Executable code
   xtmp(:) = swiftest_plA%xb(:,1)
   vtmp(:) = swiftest_plA%vb(:,1)
   do i = 1, npl
      swiftest_plA%xh(:,i) = swiftest_plA%xb(:,i) - xtmp(:)
      swiftest_plA%vh(:,i) = swiftest_plA%vb(:,i) - vtmp(:)
   end do

   return

end subroutine coord_b2h