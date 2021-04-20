subroutine coord_vh2vb(npl, swiftest_plA, msys)
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, planet velocities only
   !!  
   !! Adapted from David E. Kaufmann Swifter routine coord_vh2vb.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_vh2b.f
   use swiftest
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => coord_vh2vb
   implicit none

! arguments
   integer(I4B), intent(in)     :: npl
   real(DP), intent(out)        :: msys
   type(swiftest_pl), intent(inout) :: swiftest_plA

! internals
   integer(I4B)          :: i
   real(DP), dimension(NDIM) :: vtmp

! executable code

   if (any(swiftest_plA%status(2:npl) /= ACTIVE)) then
      vtmp(:) = 0.0_DP
      do i = 2, npl
         if (swiftest_plA%status(i) /= ACTIVE) cycle
         vtmp(:) = vtmp(:) + swiftest_plA%mass(i)*swiftest_plA%vh(:,i)
      end do
      msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl), swiftest_plA%status(2:npl) == ACTIVE) + swiftest_plA%Mcb_initial
      swiftest_plA%vb(:,1) = -vtmp(:)  /msys
      vtmp(:) = swiftest_plA%vb(:,1)
      do i = 2, npl
         if (swiftest_plA%status(i) /= ACTIVE) cycle
         swiftest_plA%vb(:,i) = swiftest_plA%vh(:,i) + vtmp(:)
      end do
   else
      vtmp(:) = matmul(swiftest_plA%vh(:,2:npl), swiftest_plA%mass(2:npl))
      msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl)) + swiftest_plA%Mcb_initial
      swiftest_plA%vb(:,1) = -vtmp(:) / msys    
      do i = 1, NDIM
         swiftest_plA%vb(i,2:npl) = swiftest_plA%vh(i,2:npl) + swiftest_plA%vb(i,1)
      end do 
   end if

   return

end subroutine coord_vh2vb