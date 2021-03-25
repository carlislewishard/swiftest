subroutine coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)
   !! Author: David Minton
   !! 
   !! Convert from barycentric to heliocentric coordinates, active test particles only
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine coord_b2h_tp.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_b2h_tp.f90
! Modules
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => coord_b2h_tp
   implicit none

! Arguments
   integer(I4B), intent(in)  :: ntp
   type(swiftest_tp), intent(inout) :: swiftest_tpA
   type(swiftest_pl), intent(inout) :: swiftest_plA

! Internals
   integer(I4B)          :: i
   real(DP), dimension(NDIM) :: xtmp, vtmp

! Executable code
   xtmp(:) = swiftest_plA%xb(:,1)
   vtmp(:) = swiftest_plA%vb(:,1)
   do i = 1, ntp
      swiftest_tpA%xh(:,i) = swiftest_tpA%xb(:,i) - xtmp(:)
      swiftest_tpA%vh(:,i) = swiftest_tpA%vb(:,i) - vtmp(:)
   end do

   return

end subroutine coord_b2h_tp