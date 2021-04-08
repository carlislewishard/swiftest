subroutine symba_rearray_pl(npl, symba_plA, nmergeadd, mergeadd_list, discard_plA)
   !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Clean up tp and pl arrays to remove discarded bodies and add new bodies
! modules
   use swiftest
   use module_swiftestalloc 
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_rearray_pl
   implicit none

! arguments
   integer(I4B), intent(inout)             :: npl,  nmergeadd 
   type(symba_pl), intent(inout)           :: symba_plA
   type(swiftest_pl), intent(inout)        :: discard_plA
   type(symba_merger), intent(inout)       :: mergeadd_list 
! internals
   integer(I4B)                           :: i, nkpl, nktp, ntot, nsppl
   integer(I4B)                           :: ip, nchild, j
   logical, dimension(:), allocatable     :: discard_l_pl 
   real(DP)                               :: msys_old, msys_new

! executable code

   nkpl = npl
   call util_resize_pl(symba_plA, npl +nmergeadd, npl)
   npl = nkpl + nmergeadd

   !add merge products to the end of the planet list
   symba_plA%helio%swiftest%status(nkpl+1:npl) = ACTIVE
   symba_plA%helio%swiftest%name(nkpl+1:npl)   = mergeadd_list%name(1:nmergeadd)
   symba_plA%helio%swiftest%mass(nkpl+1:npl)   = mergeadd_list%mass(1:nmergeadd)
   symba_plA%helio%swiftest%radius(nkpl+1:npl) = mergeadd_list%radius(1:nmergeadd)
   do i = 1, NDIM
      symba_plA%helio%swiftest%xh(i,nkpl+1:npl)  = mergeadd_list%xh(i,1:nmergeadd)
      symba_plA%helio%swiftest%vh(i,nkpl+1:npl)  = mergeadd_list%vh(i,1:nmergeadd)
      symba_plA%helio%swiftest%ip(i,nkpl+1:npl)  = mergeadd_list%ip(i,1:nmergeadd)
      symba_plA%helio%swiftest%rot(i,nkpl+1:npl) = mergeadd_list%rot(i,1:nmergeadd)
   end do

   call util_hills(npl, symba_plA%helio%swiftest)

   npl = count(symba_plA%helio%swiftest%status(1:npl) /= INACTIVE)
   ntot= size(symba_plA%helio%swiftest%status(:))
   if (ntot > npl) then
      symba_plA%helio%swiftest%status(npl+1:ntot) = INACTIVE
      symba_plA%helio%swiftest%mass(npl+1:ntot) = 0.0_DP
      symba_plA%helio%swiftest%name(npl+1:ntot) = -1
   end if

   symba_plA%helio%swiftest%nbody = npl

   !call symba_reorder_pl(npl, symba_plA)


end subroutine symba_rearray_pl
