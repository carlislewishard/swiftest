subroutine coord_h2b(npl, swiftest_plA, msys)
   !! author: David A. Minton
   !!
   !! Convert from heliocentric to barycentric coordinates, planets only
   !!  
   !! Adapted from David E. Kaufmann Swifter routine coord_h2b.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine coord_h2b.f
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => coord_h2b
   implicit none

! arguments
   integer(I4B), intent(in)  :: npl
   real(DP), intent(out)   :: msys
   type(swiftest_pl),intent(inout) :: swiftest_plA

! internals
   integer(I4B)          :: i
   logical, dimension(npl) :: lstatus

! executable code
   associate(vbcb => swiftest_plA%vb(:,1), xbcb => swiftest_plA%xb(:,1), &
             vb => swiftest_plA%vb,        vh   => swiftest_plA%vh, &
             xb => swiftest_plA%xb,        xh   => swiftest_plA%xh, &
             mass => swiftest_plA%mass, status => swiftest_plA%status, &
             dMcb => swiftest_plA%dMcb, Mcb_initial => swiftest_plA%Mcb_initial)

      lstatus(2:npl) = status(2:npl) == ACTIVE

      xbcb(:) = 0.0_DP
      vbcb(:) = 0.0_DP
      do i = 2,npl
         if (.not.lstatus(i)) cycle
         xbcb(:) = xbcb(:) + mass(i) * xh(:,i)
         vbcb(:) = vbcb(:) + mass(i) * vh(:,i)
      end do

      msys = dMcb + sum(mass(2:npl), lstatus(2:npl)) + Mcb_initial

      xbcb(:) = -xbcb(:) / msys                      
      vbcb(:) = -vbcb(:) / msys    

      do i = 2,npl 
         xb(:,i) = xh(:,i) + xbcb(:)
         vb(:,i) = vh(:,i) + vbcb(:)
      end do

      status(1) = ACTIVE

   end associate

   return

end subroutine coord_h2b