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

! executable code

   swiftest_plA%vb(:,1) = -matmul(swiftest_plA%vb(:,2:npl), swiftest_plA%mass(2:npl)) / swiftest_plA%mass(1)
   do i = 1, NDIM
      swiftest_plA%vh(i,2:npl) = swiftest_plA%vb(i,2:npl) - swiftest_plA%vb(i,1) 
   end do

   return

end subroutine coord_vb2vh