subroutine coord_vb2vh(npl, swiftest_plA)
   !! author: David A. Minton
   !!
   !! Convert from barycentric to heliocentric coordinates, planet velocities only
   !!  
   !! Adapted from David E. Kaufmann Swifter routine coord_vb2vh.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_vh2h.f
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => coord_vb2vh
   implicit none

! arguments
   integer(I4B), intent(in)  :: npl
   type(swiftest_pl), intent(inout) :: swiftest_plA

! internals
   integer(I4B)          :: i
   !real(DP), dimension(NDIM) :: vtmp

! executable code

   associate(vbcb => swiftest_plA%vb(:,1), vb => swiftest_plA%vb, mass => swiftest_plA%mass, vh => swiftest_plA%vh)
      vbcb(:) = 0.0_DP
      do i = 2, npl
         vbcb(:) = vbcb(:) - mass(i) * vb(:,i)
      end do
      vbcb(:) = vbcb(:) / mass(1)
      do i = 2, npl
         vh(:,i) = vb(:,i) - vbcb(:)
      end do
   end associate

   return

end subroutine coord_vb2vh