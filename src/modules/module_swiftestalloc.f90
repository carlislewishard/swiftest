module module_swiftestalloc
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Module containing subroutines that allocate and initialize the Swiftest data structures
   !!
   use swiftest_globals
   use swiftest_data_structures

   contains 

   subroutine helio_pl_allocate(helio_plA, npl)
      use module_helio
      implicit none

      integer(I4B), intent(in)            :: npl
      type(helio_pl), intent(inout)        :: helio_plA
      integer(I4B)                            :: n

      if (npl <= 0) then
         n = 1
      else
         n = npl
      end if

      allocate(helio_plA%ah(NDIM, n))
      allocate(helio_plA%ahi(NDIM, n))
      helio_plA%ah(:, :)= 0.0_DP
      helio_plA%ahi(:, :) = 0.0_DP
      call helio_plA%swiftest%alloc(npl)
      return
   end subroutine helio_pl_allocate


   subroutine symba_pl_allocate(symba_plA, npl)
      use module_symba
      implicit none

      integer(I4B), intent(in)            :: npl
      type(symba_pl), intent(inout)        :: symba_plA
      integer(I4B)                            :: n

      if (npl <= 0) then
         n = 1
      else
         n = npl
      end if

      allocate(symba_plA%lmerged(n))
      allocate(symba_plA%nplenc(n))
      allocate(symba_plA%ntpenc(n))
      allocate(symba_plA%levelg(n))
      allocate(symba_plA%levelm(n))
      allocate(symba_plA%nchild(n))
      allocate(symba_plA%isperi(n))
      allocate(symba_plA%peri(n))
      allocate(symba_plA%atp(n))
      allocate(symba_plA%index_parent(n))
      allocate(symba_plA%index_child(NCHILDMAX, n))

      symba_plA%lmerged(:) = .false.
      symba_plA%nplenc(:) = 0
      symba_plA%ntpenc(:) = 0
      symba_plA%levelg(:) = 0
      symba_plA%levelm(:) = 0
      symba_plA%nchild(:) = 0
      symba_plA%isperi(:) = 0
      symba_plA%peri(:) = 0.0_DP
      symba_plA%atp(:) = 0.0_DP
      symba_plA%index_parent(:) = 1
      symba_plA%index_child(:, :) = 1
      call helio_pl_allocate(symba_plA%helio,npl)
      return
   end subroutine symba_pl_allocate

   subroutine symba_plplenc_allocate(plplenc_list, nplplenc)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: nplplenc
      type(symba_plplenc), intent(inout)        :: plplenc_list
      integer(I4B)                            :: n

      if (nplplenc <= 0) then
         n = 1
      else
         n = nplplenc
      end if

      allocate(plplenc_list%lvdotr(n))
      allocate(plplenc_list%status(n))
      allocate(plplenc_list%level(n))
      allocate(plplenc_list%index1(n))
      allocate(plplenc_list%index2(n))
      allocate(plplenc_list%enc_child(n))
      allocate(plplenc_list%enc_parent(n))

      plplenc_list%lvdotr(:) = .false.
      plplenc_list%status(:) = 0
      plplenc_list%level(:) = 1
      plplenc_list%index1(:) = 1
      plplenc_list%index2(:) = 1
      plplenc_list%enc_child(:) = 1
      plplenc_list%enc_parent(:) = 1
      return
   end subroutine symba_plplenc_allocate

   subroutine symba_merger_allocate(merger_list, nmerger)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: nmerger
      type(symba_merger), intent(inout)       :: merger_list
      integer(I4B)                            :: n

      if (nmerger <= 0) then
         n = 0
      else
         n = nmerger
      end if

      allocate(merger_list%name(n))
      allocate(merger_list%index_ps(n))
      allocate(merger_list%status(n))
      allocate(merger_list%ncomp(n))
      allocate(merger_list%nadded(n))
      allocate(merger_list%xh(NDIM, n))
      allocate(merger_list%vh(NDIM, n))
      allocate(merger_list%mass(n))
      allocate(merger_list%radius(n))

      merger_list%name(:) = 0
      merger_list%index_ps(:) = 1
      merger_list%status(:) = 0
      merger_list%ncomp(:) = 0
      merger_list%nadded(:) = 0
      merger_list%xh(:, :) = 0.0_DP
      merger_list%vh(:, :) = 0.0_DP
      merger_list%mass(:) = 0.0_DP
      merger_list%radius(:) = 0.0_DP

      return
   end subroutine symba_merger_allocate

   subroutine helio_tp_allocate(helio_tpA, ntp)
      use module_helio
      implicit none

      integer(I4B), intent(in)            :: ntp
      type(helio_tp), intent(inout)       :: helio_tpA
      integer(I4B)                        :: n

      if (ntp <= 0) then
         n = 1
      else
         n = ntp
      end if

      allocate(helio_tpA%ah(NDIM, n))
      allocate(helio_tpA%ahi(NDIM, n))

      helio_tpA%ah(:, :) = 0.0_DP
      helio_tpA%ahi(:, :) = 0.0_DP
      call helio_tpA%swiftest%alloc(ntp)

      return
   end subroutine helio_tp_allocate


   subroutine symba_tp_allocate(symba_tpA, ntp)
      use module_symba
      implicit none

      integer(I4B), intent(in)            :: ntp
      type(symba_tp), intent(inout)       :: symba_tpA
      integer(I4B)                        :: n

      if (ntp <= 0) then
         n = 1
      else
         n = ntp
      end if
      allocate(symba_tpA%nplenc(n))
      allocate(symba_tpA%levelg(n))
      allocate(symba_tpA%levelm(n))

      symba_tpA%nplenc(:) = 0
      symba_tpA%levelg(:) = 0
      symba_tpA%levelm(:) = 0
      call helio_tp_allocate(symba_tpA%helio,ntp)
      return
   end subroutine symba_tp_allocate

   subroutine symba_pltpenc_allocate(pltpenc_list, npltpenc)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: npltpenc
      type(symba_pltpenc), intent(inout)        :: pltpenc_list
      integer(I4B)                            :: n

      if (npltpenc <= 0)  then
         n = 1
      else
         n = npltpenc
      end if

      allocate(pltpenc_list%lvdotr(n))
      allocate(pltpenc_list%status(n))
      allocate(pltpenc_list%level(n))
      allocate(pltpenc_list%indexpl(n))
      allocate(pltpenc_list%indextp(n))

      pltpenc_list%lvdotr(:) = .false.
      pltpenc_list%status(:) = 0
      pltpenc_list%level(:) = 0
      pltpenc_list%indexpl(:) = 1
      pltpenc_list%indextp(:) = 1
      return
   end subroutine symba_pltpenc_allocate

!___________________________


   subroutine helio_pl_deallocate(helio_plA)
      use module_helio
      implicit none

      type(helio_pl), intent(inout)        :: helio_plA

      if (allocated(helio_plA%ah)) deallocate(helio_plA%ah)
      if (allocated(helio_plA%ahi)) deallocate(helio_plA%ahi)
      call helio_plA%swiftest%dealloc()
      return
   end subroutine helio_pl_deallocate


   subroutine symba_pl_deallocate(symba_plA)
      use module_symba
      implicit none

      type(symba_pl), intent(inout)        :: symba_plA

      if (allocated(symba_plA%lmerged)) deallocate(symba_plA%lmerged)
      if (allocated(symba_plA%nplenc)) deallocate(symba_plA%nplenc)
      if (allocated(symba_plA%ntpenc)) deallocate(symba_plA%ntpenc)
      if (allocated(symba_plA%levelg)) deallocate(symba_plA%levelg)
      if (allocated(symba_plA%levelm)) deallocate(symba_plA%levelm)
      if (allocated(symba_plA%nchild)) deallocate(symba_plA%nchild)
      if (allocated(symba_plA%isperi)) deallocate(symba_plA%isperi)
      if (allocated(symba_plA%peri)) deallocate(symba_plA%peri)
      if (allocated(symba_plA%atp)) deallocate(symba_plA%atp)
      if (allocated(symba_plA%index_parent)) deallocate(symba_plA%index_parent)
      if (allocated(symba_plA%index_child)) deallocate(symba_plA%index_child)
      call helio_pl_deallocate(symba_plA%helio)
      return
   end subroutine symba_pl_deallocate

   subroutine symba_plplenc_deallocate(plplenc_list)
      use module_symba
      implicit none

      type(symba_plplenc), intent(inout)        :: plplenc_list

      if (allocated(plplenc_list%lvdotr)) deallocate(plplenc_list%lvdotr)
      if (allocated(plplenc_list%status)) deallocate(plplenc_list%status)
      if (allocated(plplenc_list%level)) deallocate(plplenc_list%level)
      if (allocated(plplenc_list%index1)) deallocate(plplenc_list%index1)
      if (allocated(plplenc_list%index2)) deallocate(plplenc_list%index2)
      if (allocated(plplenc_list%enc_child)) deallocate(plplenc_list%enc_child)
      if (allocated(plplenc_list%enc_parent)) deallocate(plplenc_list%enc_parent)
      return
   end subroutine symba_plplenc_deallocate

   subroutine symba_merger_deallocate(mergeadd_list)
      use module_symba
      implicit none

      type(symba_merger), intent(inout)        :: mergeadd_list

      if (allocated(mergeadd_list%name)) deallocate(mergeadd_list%name)
      if (allocated(mergeadd_list%index_ps)) deallocate(mergeadd_list%index_ps)
      if (allocated(mergeadd_list%status)) deallocate(mergeadd_list%status)
      if (allocated(mergeadd_list%ncomp)) deallocate(mergeadd_list%ncomp)
      if (allocated(mergeadd_list%nadded)) deallocate(mergeadd_list%nadded)
      if (allocated(mergeadd_list%xh)) deallocate(mergeadd_list%xh)
      if (allocated(mergeadd_list%vh)) deallocate(mergeadd_list%vh)
      if (allocated(mergeadd_list%mass)) deallocate(mergeadd_list%mass)
      if (allocated(mergeadd_list%radius)) deallocate(mergeadd_list%radius)
      return
   end subroutine symba_merger_deallocate

   subroutine helio_tp_deallocate(helio_tpA)
      use module_helio
      implicit none

      type(helio_tp), intent(inout)        :: helio_tpA

      if (allocated(helio_tpA%ah)) deallocate(helio_tpA%ah)
      if (allocated(helio_tpA%ahi)) deallocate(helio_tpA%ahi)
      call helio_tpA%swiftest%dealloc()
      return
   end subroutine helio_tp_deallocate


   subroutine symba_tp_deallocate(symba_tpA)
      use module_symba
      implicit none

      type(symba_tp), intent(inout)        :: symba_tpA

      if (allocated(symba_tpA%nplenc)) deallocate(symba_tpA%nplenc)
      if (allocated(symba_tpA%levelg)) deallocate(symba_tpA%levelg)
      if (allocated(symba_tpA%levelm)) deallocate(symba_tpA%levelm)

      return
   end subroutine symba_tp_deallocate

   subroutine symba_pltpenc_deallocate(pltpenc_list)
      use module_symba
      implicit none

      type(symba_pltpenc), intent(inout)        :: pltpenc_list

      if (allocated(pltpenc_list%lvdotr)) deallocate(pltpenc_list%lvdotr)
      if (allocated(pltpenc_list%status)) deallocate(pltpenc_list%status)
      if (allocated(pltpenc_list%level)) deallocate(pltpenc_list%level)
      if (allocated(pltpenc_list%indexpl)) deallocate(pltpenc_list%indexpl)
      if (allocated(pltpenc_list%indextp)) deallocate(pltpenc_list%indextp)
      return
   end subroutine symba_pltpenc_deallocate

   subroutine symba_plplenc_copy(plplenc_list_in, plplenc_list_out, n)
      use module_symba
      implicit none

      type(symba_plplenc), intent(in)        :: plplenc_list_in
      type(symba_plplenc), intent(inout)     :: plplenc_list_out
      integer(I4B), intent(in)               :: n

      plplenc_list_out%lvdotr(1:n) = plplenc_list_in%lvdotr(1:n)
      plplenc_list_out%status(1:n) = plplenc_list_in%status(1:n)
      plplenc_list_out%level(1:n) = plplenc_list_in%level(1:n)
      plplenc_list_out%index1(1:n) = plplenc_list_in%index1(1:n)
      plplenc_list_out%index2(1:n) = plplenc_list_in%index2(1:n)
      plplenc_list_out%enc_child(1:n) = plplenc_list_in%enc_child(1:n)
      plplenc_list_out%enc_parent(1:n) = plplenc_list_in%enc_parent(1:n)
      return
   end subroutine symba_plplenc_copy

   subroutine symba_pltpenc_copy(pltpenc_list_in, pltpenc_list_out, n)
      use module_symba
      implicit none

      type(symba_pltpenc), intent(in)        :: pltpenc_list_in
      type(symba_pltpenc), intent(inout)     :: pltpenc_list_out
      integer(I4B), intent(in)               :: n

      pltpenc_list_out%lvdotr(1:n) = pltpenc_list_in%lvdotr(1:n)
      pltpenc_list_out%status(1:n) = pltpenc_list_in%status(1:n)
      pltpenc_list_out%level(1:n) = pltpenc_list_in%level(1:n)
      pltpenc_list_out%indexpl(1:n) = pltpenc_list_in%indexpl(1:n)
      pltpenc_list_out%indextp(1:n) = pltpenc_list_in%indextp(1:n)

      return
   end subroutine symba_pltpenc_copy

   subroutine symba_merger_copy(merger_list_in, merger_list_out, n)
      use module_symba
      implicit none

      type(symba_merger), intent(in)        :: merger_list_in
      type(symba_merger), intent(inout)     :: merger_list_out
      integer(I4B), intent(in)              :: n

      merger_list_out%name(1:n) = merger_list_in%name(1:n)
      merger_list_out%index_ps(1:n) = merger_list_in%index_ps(1:n)
      merger_list_out%status(1:n) = merger_list_in%status(1:n)
      merger_list_out%ncomp(1:n) = merger_list_in%ncomp(1:n)
      merger_list_out%nadded(1:n) = merger_list_in%nadded(1:n)
      merger_list_out%xh(:, 1:n) = merger_list_in%xh(:, 1:n)
      merger_list_out%vh(:, 1:n) = merger_list_in%vh(:, 1:n)
      merger_list_out%mass(1:n) = merger_list_in%mass(1:n)
      merger_list_out%radius(1:n) = merger_list_in%radius(1:n)
      return
   end subroutine symba_merger_copy

   subroutine symba_plplenc_size_check(plplenc_list, n)
      !! author: David A. Minton
      !!
      !! Checks the current size of the plplenc_list against the required size and extends 
      !! it by a factor of 2 if it is too small 
      use module_symba
      implicit none

      type(symba_plplenc), intent(inout)    :: plplenc_list
      integer(I4B), intent(in)              :: n 
      type(symba_plplenc)                   :: plplenc_temp

      if (n <= size(plplenc_list%status)) return

      call symba_plplenc_allocate(plplenc_temp, n)
      call symba_plplenc_copy(plplenc_list, plplenc_temp,n - 1)
      call symba_plplenc_deallocate(plplenc_list)
      call symba_plplenc_allocate(plplenc_list, 2 * n)
      call symba_plplenc_copy(plplenc_temp, plplenc_list, n - 1)
      call symba_plplenc_deallocate(plplenc_temp)

      return

   end subroutine symba_plplenc_size_check

   subroutine symba_pltpenc_size_check(pltpenc_list, n)
      !! author: David A. Minton
      !!
      !! Checks the current size of the pltpenc_list against the required size and extends 
      !! it by a factor of 2 if it is too small 
      use module_symba
      implicit none

      type(symba_pltpenc), intent(inout)    :: pltpenc_list
      integer(I4B), intent(in)              :: n 
      type(symba_pltpenc)                   :: pltpenc_temp

      if (n <= size(pltpenc_list%status)) return

      call symba_pltpenc_allocate(pltpenc_temp, n)
      call symba_pltpenc_copy(pltpenc_list, pltpenc_temp,n - 1)
      call symba_pltpenc_deallocate(pltpenc_list)
      call symba_pltpenc_allocate(pltpenc_list, 2 * n)
      call symba_pltpenc_copy(pltpenc_temp, pltpenc_list, n - 1)
      call symba_pltpenc_deallocate(pltpenc_temp)

      return

   end subroutine symba_pltpenc_size_check

   subroutine symba_merger_size_check(merger_list, n)
      !! author: David A. Minton
      !!
      !! Checks the current size of the merger_list against the required size and extends 
      !! it by a factor of 2 if it is too small 
      use module_symba
      implicit none

      type(symba_merger), intent(inout)   :: merger_list
      integer(I4B), intent(in)            :: n 
      type(symba_merger)                  :: merger_temp

      if (n <= size(merger_list%status)) return

      call symba_merger_allocate(merger_temp, n)
      call symba_merger_copy(merger_list, merger_temp,n - 1)
      call symba_merger_deallocate(merger_list)
      call symba_merger_allocate(merger_list, 2 * n)
      call symba_merger_copy(merger_temp, merger_list, n - 1)
      call symba_merger_deallocate(merger_temp)

      return

   end subroutine symba_merger_size_check


end module module_swiftestalloc




