subroutine symba_rearray(npl, nplm, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA,&
   discard_tpA, ldiscard_pl, ldiscard_tp, mtiny)
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
   integer(I4B), intent(inout)             :: npl, nplm, ntp, nsppl, nsptp
   integer(I4B), intent(in)                :: nmergeadd 
   type(symba_pl), intent(inout)           :: symba_plA
   type(symba_tp), intent(inout)           :: symba_tpA
   type(symba_tp), intent(inout)           :: discard_tpA
   type(symba_pl), intent(inout)           :: discard_plA
   type(symba_merger), intent(inout)       :: mergeadd_list 
   logical, intent(in)                     :: ldiscard_pl, ldiscard_tp 
   real(DP), intent(in)                    :: mtiny

   ! Internals
   integer(I4B)                           :: i, j, nkpl, nktp, ntot, dlo
   logical, dimension(:), allocatable     :: discard_l_pl 
   logical, dimension(ntp)                :: discard_l_tp
   logical                                :: lescape

   if (any(symba_plA%helio%swiftest%status(1:npl) /= ACTIVE)) then 

      ! Deal with the central body/system discards if there are any
      do i = 1, npl
         if ((symba_plA%helio%swiftest%status(i) == DISCARDED_RMIN) .or. (symba_plA%helio%swiftest%status(i) == DISCARDED_PERI)) then
            lescape = .false.
         else if ((symba_plA%helio%swiftest%status(i) == DISCARDED_RMAX) .or. (symba_plA%helio%swiftest%status(i) == DISCARDED_RMAXU)) then
            lescape = .true.
         else 
            cycle
         end if
         call symba_discard_conserve_mtm(symba_plA%helio%swiftest, i, lescape)
      end do

      allocate(discard_l_pl(npl))
      nkpl = count(symba_plA%helio%swiftest%status(:) == ACTIVE)
      nsppl = npl - nkpl
      dlo = 1
      if (allocated(discard_plA%helio%swiftest%id)) then ! We alredy made a discard list in this step, so we need to append to it
         nsppl = nsppl + discard_plA%helio%swiftest%nbody
         dlo = dlo + discard_plA%helio%swiftest%nbody
         call util_resize_pl(discard_plA, nsppl)
      else
         call symba_pl_allocate(discard_plA, nsppl)
      end if
      !write(*,*) 'Keeping ',nkpl,' out of ',npl,' bodies.'

      discard_l_pl(:) = (symba_plA%helio%swiftest%status(1:npl) /= ACTIVE) 
      !write(*,*) 'Discarding ',count(discard_l_pl(:)), ' out of ',npl,' bodies.'
      !do i = 1, npl
      !   if (discard_l_pl(i)) write(*,*) i,symba_plA%helio%swiftest%id(i), ' is on its way to destruction. Mass: ',symba_plA%helio%swiftest%mass(i)
      !end do

      ! Spill discarded bodies into discard list
      discard_plA%helio%swiftest%id(dlo:nsppl)   = pack(symba_plA%helio%swiftest%id(1:npl),   discard_l_pl)
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
      !write(*,*) 'Total mass of keeps:    ',sum(symba_plA%helio%swiftest%mass(1:nkpl))
   else
      nkpl = npl
   end if

   if (nmergeadd > 0) then
      call util_resize_pl(symba_plA, nkpl+nmergeadd)
      npl = nkpl + nmergeadd
      symba_plA%helio%swiftest%nbody = npl

      !write(*,*) 'Adding ',nmergeadd, ' bodies.'
      !do i = 1, nmergeadd 
      !   write(*,*) nkpl+i,mergeadd_list%id(i), ' is on its way to creation.    Mass: ', mergeadd_list%mass(i)
      !end do

      !add merge products to the end of the planet list
      symba_plA%helio%swiftest%status(nkpl+1:npl) = ACTIVE
      symba_plA%helio%swiftest%id(nkpl+1:npl)   = mergeadd_list%id(1:nmergeadd)
      symba_plA%helio%swiftest%mass(nkpl+1:npl)   = mergeadd_list%mass(1:nmergeadd)
      symba_plA%helio%swiftest%radius(nkpl+1:npl) = mergeadd_list%radius(1:nmergeadd)
      do i = 1, NDIM
         symba_plA%helio%swiftest%xb(i,nkpl+1:npl)  = mergeadd_list%xb(i,1:nmergeadd)
         symba_plA%helio%swiftest%vb(i,nkpl+1:npl)  = mergeadd_list%vb(i,1:nmergeadd)
         symba_plA%helio%swiftest%Ip(i,nkpl+1:npl)  = mergeadd_list%Ip(i,1:nmergeadd)
         symba_plA%helio%swiftest%rot(i,nkpl+1:npl) = mergeadd_list%rot(i,1:nmergeadd)
      end do
      !write(*,*) 'Total mass of discards: ',sum(symba_plA%helio%swiftest%mass(nkpl+1:npl))
   end if 

   if (ldiscard_pl) then
      ntot= size(symba_plA%helio%swiftest%status(:))
      if (ntot > npl) then
         symba_plA%helio%swiftest%status(npl+1:ntot) = INACTIVE
         symba_plA%helio%swiftest%mass(npl+1:ntot) = 0.0_DP
         symba_plA%helio%swiftest%id(npl+1:ntot) = -1
      end if

      call symba_reorder_pl(npl, symba_plA)

      nplm = count(symba_plA%helio%swiftest%mass > mtiny)

      call coord_b2h(npl, symba_plA%helio%swiftest)
      
      call util_hills(npl, symba_plA%helio%swiftest)

      CALL util_dist_index_plpl(npl, nplm, symba_plA)

      ! check for duplicate names and fix if ncessary
      do j = 1, npl
         do i = j + 1, npl
            if (symba_plA%helio%swiftest%id(i) == symba_plA%helio%swiftest%id(j)) then
               symba_plA%helio%swiftest%id(i) = maxval(symba_plA%helio%swiftest%id(:)) + 1
            end if
         end do
      end do
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
