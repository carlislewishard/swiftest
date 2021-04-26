subroutine symba_reorder_pl(npl, symba_plA)
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Rearrange symba planet arrays in order of decreasing mass
   !!   
   use swiftest
   use swiftest_globals
   use swiftest_data_structures
   use module_swiftestalloc
   USE module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_reorder_pl
   implicit none

   integer(I4B), intent(in) :: npl
   type(symba_pl), intent(inout)  :: symba_plA

   integer(I4B)                  :: i,j,narr
   integer(I4B), dimension(:), allocatable   :: sort_index
   type(symba_pl)                  :: symba_plwkspA

   call symba_pl_allocate(symba_plwkspA,npl)
   allocate(sort_index(npl))

   ! sort by mass
   call util_index(symba_plA%helio%swiftest%mass(2:npl), sort_index(2:npl))
   !write(*,*) "************ REORDER ***************"
   do i = 2, npl
      symba_plA%helio%swiftest%name(i) = symba_plwkspA%helio%swiftest%name(sort_index(npl-i+2))
      symba_plA%helio%swiftest%status(i) = symba_plwkspA%helio%swiftest%status(sort_index(npl-i+2))
      symba_plA%helio%swiftest%mass(i) = symba_plwkspA%helio%swiftest%mass(sort_index(npl-i+2))
      symba_plA%helio%swiftest%radius(i) = symba_plwkspA%helio%swiftest%radius(sort_index(npl-i+2))
      symba_plA%helio%swiftest%xh(:,i) = symba_plwkspA%helio%swiftest%xh(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%vh(:,i) = symba_plwkspA%helio%swiftest%vh(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%xb(:,i) = symba_plwkspA%helio%swiftest%xb(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%vb(:,i) = symba_plwkspA%helio%swiftest%vb(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%rot(:,i) = symba_plwkspA%helio%swiftest%rot(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%Ip(:,i) = symba_plwkspA%helio%swiftest%Ip(:,sort_index(npl-i+2))
      symba_plA%helio%swiftest%rhill(i) = symba_plwkspA%helio%swiftest%rhill(sort_index(npl-i+2))
   end do
   
   narr = size(symba_plA%helio%swiftest%status(:))
   symba_plA%helio%swiftest%name(npl+1:narr) = 0
   symba_plA%helio%swiftest%status(npl+1:narr) = INACTIVE
   call symba_pl_deallocate(symba_plwkspA)
   deallocate(sort_index)

   return

end subroutine symba_reorder_pl