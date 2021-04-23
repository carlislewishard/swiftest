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
   logical, dimension(npl) :: lstatus

! executable code

   associate(vbcb => swiftest_plA%vb(:,1), &
             vb   => swiftest_plA%vb,   vh => swiftest_plA%vh, &
             mass => swiftest_plA%mass, status => swiftest_plA%status, &
             dMcb => swiftest_plA%dMcb, Mcb_initial => swiftest_plA%Mcb_initial)

      lstatus(2:npl) = status(2:npl) == ACTIVE

      vbcb(:) = 0.0_DP
      do i = 2,npl
         if (.not.lstatus(i)) cycle
         vbcb(:) = vbcb(:) + mass(i) * vh(:,i)
      end do

      msys = dMcb + sum(mass(2:npl), lstatus(2:npl)) + Mcb_initial
      vbcb(:) = -vbcb(:) / msys

      do i = 2, npl
         if (.not.lstatus(i)) cycle
         vb(:,i) = vh(:,i) + vbcb(:)
      end do

   end associate

   return

end subroutine coord_vh2vb