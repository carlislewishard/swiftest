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
   logical, dimension(npl) :: lstatus

! Executable code

   associate(vbcb => swiftest_plA%vb(:,1), xbcb => swiftest_plA%xb(:,1), &
      vb => swiftest_plA%vb,        vh   => swiftest_plA%vh, &
      xb => swiftest_plA%xb,        xh   => swiftest_plA%xh)

      lstatus(2:npl) = status(2:npl) == ACTIVE

      do i = 1, npl
         if (.not.lstatus(i)) cycle
         xh(:,i) = xb(:,i) - xbcb(:)
         vh(:,i) = vb(:,i) - vbcb(:)
      end do

   end associate

   return

end subroutine coord_b2h