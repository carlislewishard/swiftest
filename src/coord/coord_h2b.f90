subroutine coord_h2b(npl, swiftest_plA, msys)
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, planets only
   !!  
   !! Adapted from David E. Kaufmann Swifter routine coord_h2b.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_h2b.f
   use swiftest
   use module_interfaces, except_this_one => coord_h2b
   implicit none

! arguments
   integer(I4B), intent(in)  :: npl
   real(DP), intent(out)   :: msys
   type(swiftest_pl),intent(inout) :: swiftest_plA

! internals
   integer(I4B)          :: i
   real(DP), dimension(NDIM) :: xtmp, vtmp

! executable code
   if (any(swiftest_plA%status(2:npl) /= ACTIVE)) then
      xtmp(:) = 0.0_DP
      vtmp(:) = 0.0_DP
      do i = 2, npl
         if (swiftest_plA%status(i) /= ACTIVE) cycle
         xtmp(:) = xtmp(:) + swiftest_plA%mass(i)*swiftest_plA%xh(:,i)
         vtmp(:) = vtmp(:) + swiftest_plA%mass(i)*swiftest_plA%vh(:,i)
      end do
      msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl), swiftest_plA%status(2:npl) == ACTIVE) + swiftest_plA%Mcb_initial
      swiftest_plA%xb(:,1) = -xtmp(:) / msys                      
      swiftest_plA%vb(:,1) = -vtmp(:) / msys                      
      xtmp(:) = swiftest_plA%xb(:,1)
      vtmp(:) = swiftest_plA%vb(:,1)
      do i = 2, npl
         if (swiftest_plA%status(i) /= ACTIVE) cycle
         swiftest_plA%xb(:,i) = swiftest_plA%xh(:,i) + xtmp(:)
         swiftest_plA%vb(:,i) = swiftest_plA%vh(:,i) + vtmp(:)
      end do
   else
      xtmp(:) = matmul(swiftest_plA%xh(:,2:npl), swiftest_plA%mass(2:npl))
      vtmp(:) = matmul(swiftest_plA%vh(:,2:npl), swiftest_plA%mass(2:npl))
      msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl)) + swiftest_plA%Mcb_initial
      swiftest_plA%xb(:,1) = -xtmp(:) / msys                      
      swiftest_plA%vb(:,1) = -vtmp(:) / msys    
      do i = 1, NDIM
         swiftest_plA%xb(i,2:npl) = swiftest_plA%xh(i,2:npl) + swiftest_plA%xb(i,1)
         swiftest_plA%vb(i,2:npl) = swiftest_plA%vh(i,2:npl) + swiftest_plA%vb(i,1)
      end do

   end if

   return

end subroutine coord_h2b