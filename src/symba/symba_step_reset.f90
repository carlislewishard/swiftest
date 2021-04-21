subroutine symba_step_reset(npl, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list)
   !! author: David A. Minton
   !!
   !! Resets pl, tp,and encounter structures at the start of a new step
   !!   
   use swiftest
   use swiftest_globals
   use swiftest_data_structures
   USE module_symba
   use module_swiftestalloc
   use module_interfaces, EXCEPT_THIS_ONE => symba_step_reset
   implicit none

   integer(I4B), intent(in) :: npl
   type(symba_pl), intent(inout)  :: symba_plA
   type(symba_tp), intent(inout)  :: symba_tpA
   type(symba_plplenc), intent(inout)   :: plplenc_list
   type(symba_pltpenc), intent(inout)   :: pltpenc_list
   type(symba_merger), intent(inout)    :: mergeadd_list, mergesub_list

   integer(I4B)                  :: i

   symba_plA%lcollision(:) = .false.
   symba_plA%kin(:)%parent = (/ (i, i=1, size(symba_plA%kin(:))) /)
   symba_plA%kin(:)%nchild = 0
   do i = 1, size(symba_plA%kin(:))
      if (allocated(symba_plA%kin(i)%child)) deallocate(symba_plA%kin(i)%child)
   end do
   symba_plA%nplenc(:) = 0
   symba_plA%ntpenc(:) = 0
   symba_plA%levelg(:) = 0
   symba_plA%levelm(:) = 0
   symba_tpA%nplenc(:) = 0 
   symba_tpA%levelg(:) = 0
   symba_tpA%levelm(:) = 0
   symba_plA%l_plpl_encounter = .false.

   !************************
   ! Placeholder for when the data structures get re-done
   plplenc_list%nplplenc = 0
   pltpenc_list%npltpenc = 0
   !************************
   if (.not. allocated(plplenc_list%status)) call symba_plplenc_allocate(plplenc_list, 1)
   if (.not. allocated(pltpenc_list%status)) call symba_pltpenc_allocate(pltpenc_list, 1)
   if (.not. allocated(mergeadd_list%status)) call symba_merger_allocate(mergeadd_list, 1)
   if (.not. allocated(mergesub_list%status)) call symba_merger_allocate(mergesub_list, 1)
   plplenc_list%status(:) = 0
   pltpenc_list%status(:) = 0
   mergeadd_list%status(:) = 0
   mergesub_list%status(:) = 0

   return
end subroutine symba_step_reset