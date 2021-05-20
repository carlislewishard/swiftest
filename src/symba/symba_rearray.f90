subroutine symba_rearray(t, npl, nplm, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA,&
   discard_tpA, ldiscard_pl, ldiscard_tp, mtiny, param, discard_l_pl, discard_stat_list)
   !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Clean up tp and pl arrays to remove discarded bodies and add new bodies
! modules
   use swiftest
   use module_swiftestalloc 
   use module_helio
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_rearray
   implicit none

   ! Arguments
   real(DP),                    intent(in)    :: t 
   integer(I4B),                intent(inout) :: npl, nplm, ntp, nsppl, nsptp
   integer(I4B),                intent(in)    :: nmergeadd 
   type(symba_pl),              intent(inout) :: symba_plA
   type(symba_tp),              intent(inout) :: symba_tpA
   type(symba_tp),              intent(inout) :: discard_tpA
   type(symba_pl),              intent(inout) :: discard_plA
   type(symba_merger),          intent(inout) :: mergeadd_list 
   logical,                     intent(in)    :: ldiscard_pl, ldiscard_tp 
   real(DP),                    intent(in)    :: mtiny
   type(user_input_parameters), intent(in)    :: param
   logical, dimension(:), allocatable, intent(inout) :: discard_l_pl
   integer(I4B), dimension(:), allocatable, intent(inout) :: discard_stat_list

   ! Internals
   integer(I4B)                            :: i, j, nkpl, nktp, ntot, dlo, ndiscard, dstat
   logical, dimension(:), allocatable      :: add_l_pl
   logical, dimension(ntp)                 :: discard_l_tp
   logical                                 :: lescape
   integer(I4B), dimension(:), allocatable :: discard_index_list
   real(DP)                                :: msys

   ! First resolve the central body/system discards if there are any and mark them for discard
   if (any(discard_l_pl(1:npl))) then 
      call coord_h2b(npl, symba_plA%helio%swiftest, msys)
      ndiscard = count(discard_l_pl(1:npl))
      allocate(discard_index_list(ndiscard))
      discard_index_list(:) = pack([(i, i = 1, npl)], discard_l_pl(1:npl))
      do i = 1, ndiscard
         dstat = discard_stat_list(i)
         if ((dstat == DISCARDED_RMIN) .or. (dstat == DISCARDED_PERI)) then
            lescape = .false.
         else if ((dstat == DISCARDED_RMAX) .or. (dstat == DISCARDED_RMAXU)) then
            lescape = .true.
         else 
            cycle
         end if
         ! Resolve the discard
         call symba_discard_conserve_mtm(symba_plA%helio%swiftest, discard_index_list(i), lescape)
         ! Flip the main status flag to the discard state
      end do
      symba_plA%helio%swiftest%status(discard_index_list(:)) = discard_stat_list(:)
   end if

   ! Next, remove all the bodies marked for discard from the main symba_plA structure and pack them into the discard_plA structure
   nkpl = count(symba_plA%helio%swiftest%status(:) == ACTIVE)
   nsppl = npl - nkpl
   if (nsppl > 0) then
      dlo = 1
      if (allocated(discard_plA%helio%swiftest%id)) then ! We alredy made a discard list in this step, so we need to append to it
         nsppl = nsppl + discard_plA%helio%swiftest%nbody
         dlo = dlo + discard_plA%helio%swiftest%nbody
         call util_resize_pl(discard_plA, nsppl)
      else
         call symba_pl_allocate(discard_plA, nsppl)
      end if

      discard_l_pl(:) = (symba_plA%helio%swiftest%status(1:npl) /= ACTIVE) 

      ! Spill discarded bodies into discard list
      discard_plA%helio%swiftest%id(dlo:nsppl)     = pack(symba_plA%helio%swiftest%id(1:npl),   discard_l_pl)
      discard_plA%helio%swiftest%status(dlo:nsppl) = pack(symba_plA%helio%swiftest%status(1:npl), discard_l_pl)
      discard_plA%helio%swiftest%mass(dlo:nsppl)   = pack(symba_plA%helio%swiftest%mass(1:npl),   discard_l_pl)
      discard_plA%helio%swiftest%radius(dlo:nsppl) = pack(symba_plA%helio%swiftest%radius(1:npl), discard_l_pl)
      discard_plA%helio%swiftest%rhill(dlo:nsppl)  = pack(symba_plA%helio%swiftest%rhill(1:npl),  discard_l_pl)
      do i = 1, NDIM
         discard_plA%helio%swiftest%xh(i,dlo:nsppl)  = pack(symba_plA%helio%swiftest%xh(i,1:npl),  discard_l_pl)
         discard_plA%helio%swiftest%vh(i,dlo:nsppl)  = pack(symba_plA%helio%swiftest%vh(i,1:npl),  discard_l_pl)
         discard_plA%helio%swiftest%xb(i,dlo:nsppl)  = pack(symba_plA%helio%swiftest%xb(i,1:npl),  discard_l_pl)
         discard_plA%helio%swiftest%vb(i,dlo:nsppl)  = pack(symba_plA%helio%swiftest%vb(i,1:npl),  discard_l_pl)
         discard_plA%helio%swiftest%Ip(i,dlo:nsppl)  = pack(symba_plA%helio%swiftest%Ip(i,1:npl),  discard_l_pl)
         discard_plA%helio%swiftest%rot(i,dlo:nsppl) = pack(symba_plA%helio%swiftest%rot(i,1:npl), discard_l_pl)
      end do

      if (nkpl > 0) then
         ! Pack kept bodies down 
         symba_plA%helio%swiftest%id(1:nkpl)   = pack(symba_plA%helio%swiftest%id(1:npl),   .not. discard_l_pl)
         symba_plA%helio%swiftest%status(1:nkpl) = pack(symba_plA%helio%swiftest%status(1:npl), .not. discard_l_pl)
         symba_plA%helio%swiftest%mass(1:nkpl)   = pack(symba_plA%helio%swiftest%mass(1:npl),   .not. discard_l_pl)
         symba_plA%helio%swiftest%radius(1:nkpl) = pack(symba_plA%helio%swiftest%radius(1:npl), .not. discard_l_pl)
         symba_plA%helio%swiftest%rhill(1:nkpl)  = pack(symba_plA%helio%swiftest%rhill(1:npl),  .not. discard_l_pl)
         do i = 1, NDIM
            symba_plA%helio%swiftest%xh(i,1:nkpl)  = pack(symba_plA%helio%swiftest%xh(i,1:npl),  .not. discard_l_pl)
            symba_plA%helio%swiftest%vh(i,1:nkpl)  = pack(symba_plA%helio%swiftest%vh(i,1:npl),  .not. discard_l_pl)
            symba_plA%helio%swiftest%xb(i,1:nkpl)  = pack(symba_plA%helio%swiftest%xb(i,1:npl),  .not. discard_l_pl)
            symba_plA%helio%swiftest%vb(i,1:nkpl)  = pack(symba_plA%helio%swiftest%vb(i,1:npl),  .not. discard_l_pl)
            symba_plA%helio%swiftest%Ip(i,1:nkpl)  = pack(symba_plA%helio%swiftest%Ip(i,1:npl),  .not. discard_l_pl)
            symba_plA%helio%swiftest%rot(i,1:nkpl) = pack(symba_plA%helio%swiftest%rot(i,1:npl), .not. discard_l_pl)
         end do
      end if

      npl = nkpl
   else
      nkpl = npl
   end if

   if (nmergeadd > 0) then
      call util_resize_pl(symba_plA, nkpl+nmergeadd)
      npl = nkpl + nmergeadd
      symba_plA%helio%swiftest%nbody = npl

      !add merge products to the end of the planet list
      symba_plA%helio%swiftest%status(nkpl+1:npl) = mergeadd_list%status(1:nmergeadd)
      symba_plA%helio%swiftest%id(nkpl+1:npl)   = mergeadd_list%id(1:nmergeadd)
      symba_plA%helio%swiftest%mass(nkpl+1:npl)   = mergeadd_list%mass(1:nmergeadd)
      symba_plA%helio%swiftest%radius(nkpl+1:npl) = mergeadd_list%radius(1:nmergeadd)
      do i = 1, NDIM
         symba_plA%helio%swiftest%xb(i,nkpl+1:npl)  = mergeadd_list%xb(i,1:nmergeadd)
         symba_plA%helio%swiftest%vb(i,nkpl+1:npl)  = mergeadd_list%vb(i,1:nmergeadd)
         symba_plA%helio%swiftest%Ip(i,nkpl+1:npl)  = mergeadd_list%Ip(i,1:nmergeadd)
         symba_plA%helio%swiftest%rot(i,nkpl+1:npl) = mergeadd_list%rot(i,1:nmergeadd)
      end do
      symba_plA%helio%swiftest%info(nkpl+1:npl) = mergeadd_list%info(1:nmergeadd)
   end if 

   if (ldiscard_pl) then
      call coord_b2h(npl, symba_plA%helio%swiftest)
      ! Create the particle information and set the status flags of all new particles
      allocate(add_l_pl(npl))
      
      where ((symba_plA%helio%swiftest%status(:) == DISRUPTION) .or. &
            (symba_plA%helio%swiftest%status(:) == SUPERCATASTROPHIC) .or. &
            (symba_plA%helio%swiftest%status(:) == HIT_AND_RUN))
         symba_plA%helio%swiftest%info(:)%origin_time = t
         symba_plA%helio%swiftest%info(:)%origin_xh(1) = symba_plA%helio%swiftest%xh(1,:)
         symba_plA%helio%swiftest%info(:)%origin_xh(2) = symba_plA%helio%swiftest%xh(2,:)
         symba_plA%helio%swiftest%info(:)%origin_xh(3) = symba_plA%helio%swiftest%xh(3,:)
         symba_plA%helio%swiftest%info(:)%origin_vh(1) = symba_plA%helio%swiftest%vh(1,:)
         symba_plA%helio%swiftest%info(:)%origin_vh(2) = symba_plA%helio%swiftest%vh(2,:)
         symba_plA%helio%swiftest%info(:)%origin_vh(3) = symba_plA%helio%swiftest%vh(3,:)
         add_l_pl(:) = .true.
      elsewhere
         add_l_pl(:) = .false.
      end where

      ! check for duplicate names and fix if ncessary
      do j = 1, npl
         do i = j + 1, npl
            if (symba_plA%helio%swiftest%id(i) == symba_plA%helio%swiftest%id(j)) then
               symba_plA%helio%swiftest%id(i) = maxval(symba_plA%helio%swiftest%id(:)) + 1
            end if
         end do
      end do
      symba_plA%helio%swiftest%maxid = max(symba_plA%helio%swiftest%maxid, maxval(symba_plA%helio%swiftest%id))

      call io_write_particle_pl(symba_plA%helio%swiftest, pack([(i, i=1,npl)], add_l_pl(:)), param)

      symba_plA%helio%swiftest%status(1:npl) = ACTIVE
      call symba_reorder_pl(npl, symba_plA)
      
      call util_hills(npl, symba_plA%helio%swiftest)

      nplm = count(symba_plA%helio%swiftest%mass(1:npl) > mtiny)
      CALL util_dist_index_plpl(npl, nplm, symba_plA)

   end if

   if (ldiscard_tp) then 
      nktp = 0
      nsptp = 0  

      discard_l_tp(1:ntp) = (symba_tpA%helio%swiftest%status(1:ntp) /= ACTIVE)
      nsptp = count(discard_l_tp)
      nktp = ntp - nsptp

      dlo = 1
      if (allocated(discard_tpA%helio%swiftest%id)) then ! We alredy made a discard list in this step, so we need to append to it
         nsptp = nsptp + discard_tpA%helio%swiftest%nbody
         dlo = dlo + discard_tpA%helio%swiftest%nbody
         !call util_resize_tp(discard_plA, nsppl, discard_tplA%helio%swiftest%nbody) !TODO: Implement this. Probably best to wait until this gets OOFed
      else
         call symba_tp_allocate(discard_tpA, nsptp)
      end if

      discard_tpA%helio%swiftest%id(1:nsptp)   = pack(symba_tpA%helio%swiftest%id(1:ntp),   discard_l_tp)
      discard_tpA%helio%swiftest%status(1:nsptp) = pack(symba_tpA%helio%swiftest%status(1:ntp), discard_l_tp)
      discard_tpA%helio%swiftest%isperi(1:nsptp) = pack(symba_tpA%helio%swiftest%isperi(1:ntp), discard_l_tp)
      discard_tpA%helio%swiftest%peri(1:nsptp)   = pack(symba_tpA%helio%swiftest%peri(1:ntp),   discard_l_tp)
      discard_tpA%helio%swiftest%atp(1:nsptp)    = pack(symba_tpA%helio%swiftest%atp(1:ntp),    discard_l_tp)
      do i = 1, NDIM
         discard_tpA%helio%swiftest%xh(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%xh(i,1:ntp),   discard_l_tp)
         discard_tpA%helio%swiftest%vh(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%vh(i,1:ntp),   discard_l_tp)
         discard_tpA%helio%swiftest%xb(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%xb(i,1:ntp),   discard_l_tp)
         discard_tpA%helio%swiftest%vb(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%vb(i,1:ntp),   discard_l_tp)
      end do

      symba_tpA%helio%swiftest%id(1:nktp)   = pack(symba_tpA%helio%swiftest%id(1:ntp),   .not. discard_l_tp)
      symba_tpA%helio%swiftest%status(1:nktp) = pack(symba_tpA%helio%swiftest%status(1:ntp), .not. discard_l_tp)
      symba_tpA%helio%swiftest%isperi(1:nktp) = pack(symba_tpA%helio%swiftest%isperi(1:ntp), .not. discard_l_tp)
      symba_tpA%helio%swiftest%peri(1:nktp)   = pack(symba_tpA%helio%swiftest%peri(1:ntp),   .not. discard_l_tp)
      symba_tpA%helio%swiftest%atp(1:nktp)    = pack(symba_tpA%helio%swiftest%atp(1:ntp),    .not. discard_l_tp)
      do i = 1, NDIM
         symba_tpA%helio%swiftest%xh(i,1:nktp)   = pack(symba_tpA%helio%swiftest%xh(i,1:ntp), .not. discard_l_tp)
         symba_tpA%helio%swiftest%vh(i,1:nktp)   = pack(symba_tpA%helio%swiftest%vh(i,1:ntp), .not. discard_l_tp)
         symba_tpA%helio%swiftest%xb(i,1:nktp)   = pack(symba_tpA%helio%swiftest%xb(i,1:ntp), .not. discard_l_tp)
         symba_tpA%helio%swiftest%vb(i,1:nktp)   = pack(symba_tpA%helio%swiftest%vb(i,1:ntp), .not. discard_l_tp)
         symba_tpA%helio%ah(i,1:nktp)            = pack(symba_tpA%helio%ah(i,1:ntp),          .not. discard_l_tp)
      end do
      ntp = nktp
      symba_tpA%helio%swiftest%nbody = ntp
      nplm = count(symba_plA%helio%swiftest%mass>mtiny)
      call coord_b2h_tp(ntp, symba_tpA%helio%swiftest, symba_plA%helio%swiftest)
   end if 

end subroutine symba_rearray
