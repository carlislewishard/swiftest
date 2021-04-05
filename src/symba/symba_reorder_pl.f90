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
   logical, dimension(:), allocatable :: lactive

   narr = size(symba_plA%helio%swiftest%status(:))
   call symba_pl_allocate(symba_plwkspA,npl)
   allocate(lactive(narr))
   allocate(sort_index(npl))
   lactive(:) = symba_plA%helio%swiftest%status(:) == ACTIVE

   symba_plwkspA%helio%swiftest%name(:) = pack(symba_plA%helio%swiftest%name(:), lactive)
   symba_plwkspA%helio%swiftest%status(:) = pack(symba_plA%helio%swiftest%status(:), lactive)
   symba_plwkspA%helio%swiftest%mass(:) = pack(symba_plA%helio%swiftest%mass(:), lactive)
   symba_plwkspA%helio%swiftest%radius(:) = pack(symba_plA%helio%swiftest%radius(:), lactive)
   symba_plwkspA%helio%swiftest%rhill(:) = pack(symba_plA%helio%swiftest%rhill(:), lactive)
   do i = 1, NDIM
      symba_plwkspA%helio%swiftest%xh(i,:) = pack(symba_plA%helio%swiftest%xh(i,:), lactive)
      symba_plwkspA%helio%swiftest%vh(i,:) = pack(symba_plA%helio%swiftest%vh(i,:), lactive)
      symba_plwkspA%helio%swiftest%xb(i,:) = pack(symba_plA%helio%swiftest%xb(i,:), lactive)
      symba_plwkspA%helio%swiftest%vb(i,:) = pack(symba_plA%helio%swiftest%vb(i,:), lactive)
      symba_plwkspA%helio%swiftest%rot(i,:) = pack(symba_plA%helio%swiftest%rot(i,:), lactive)
      symba_plwkspA%helio%swiftest%Ip(i,:) = pack(symba_plA%helio%swiftest%Ip(i,:), lactive)
   end do

   ! sort by mass
   call util_index(pack(symba_plA%helio%swiftest%mass(:), lactive), sort_index)
   !write(*,*) "************ REORDER ***************"
   do i = 1, npl
      symba_plA%helio%swiftest%name(i) = symba_plwkspA%helio%swiftest%name(sort_index(npl-i+1))
      symba_plA%helio%swiftest%status(i) = symba_plwkspA%helio%swiftest%status(sort_index(npl-i+1))
      symba_plA%helio%swiftest%mass(i) = symba_plwkspA%helio%swiftest%mass(sort_index(npl-i+1))
      symba_plA%helio%swiftest%radius(i) = symba_plwkspA%helio%swiftest%radius(sort_index(npl-i+1))
      symba_plA%helio%swiftest%xh(:,i) = symba_plwkspA%helio%swiftest%xh(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%vh(:,i) = symba_plwkspA%helio%swiftest%vh(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%xb(:,i) = symba_plwkspA%helio%swiftest%xb(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%vb(:,i) = symba_plwkspA%helio%swiftest%vb(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%rot(:,i) = symba_plwkspA%helio%swiftest%rot(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%Ip(:,i) = symba_plwkspA%helio%swiftest%Ip(:,sort_index(npl-i+1))
      symba_plA%helio%swiftest%rhill(i) = symba_plwkspA%helio%swiftest%rhill(sort_index(npl-i+1))
   end do
   symba_plA%helio%swiftest%name(npl+1:narr) = 0
   symba_plA%helio%swiftest%status(npl+1:narr) = INACTIVE
   call symba_pl_deallocate(symba_plwkspA)
   deallocate(sort_index)

   return

end subroutine symba_reorder_pl