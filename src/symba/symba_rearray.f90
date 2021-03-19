subroutine symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
   discard_tpA,param, ldiscard, ldiscard_tp)
   !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Clean up tp and pl arrays to remove discarded bodies and add new bodies
! modules
   use swiftest
   use module_swiftestalloc 
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => symba_rearray
   implicit none

! arguments
   integer(I4B), intent(inout)             :: npl, ntp, nsppl, nsptp, nmergeadd 
   type(symba_pl), intent(inout)           :: symba_plA
   type(symba_tp), intent(inout)           :: symba_tpA
   type(swiftest_tp), intent(inout)        :: discard_tpA
   type(swiftest_pl), intent(inout)        :: discard_plA
   type(symba_merger), intent(inout)       :: mergeadd_list 
   type(user_input_parameters),intent(in)  :: param
   logical, intent(in)                     :: ldiscard, ldiscard_tp 

! internals
   integer(I4B)                           :: i, nkpl, nktp
   logical, dimension(npl)                :: discard_l_pl 
   logical, dimension(ntp)                :: discard_l_tp

! executable code
   if (ldiscard) then 
      nsppl = 0
      nkpl = 0
      discard_l_pl(1:npl) = (symba_plA%helio%swiftest%status(1:npl) /= ACTIVE) 
      nsppl = count(discard_l_pl)
      nkpl = npl - nsppl

      call discard_plA%alloc(nsppl)

      ! Spill discarded bodies into discard list
      discard_plA%name(1:nsppl)   = pack(symba_plA%helio%swiftest%name(1:npl), discard_l_pl)
      discard_plA%status(1:nsppl) = pack(symba_plA%helio%swiftest%status(1:npl), discard_l_pl)
      discard_plA%mass(1:nsppl)   = pack(symba_plA%helio%swiftest%mass(1:npl), discard_l_pl)
      discard_plA%radius(1:nsppl) = pack(symba_plA%helio%swiftest%radius(1:npl), discard_l_pl)
      discard_plA%rhill(1:nsppl)  = pack(symba_plA%helio%swiftest%rhill(1:npl), discard_l_pl)
      do i = 1, NDIM
         discard_plA%xh(i,1:nsppl)  = pack(symba_plA%helio%swiftest%xh(i,1:npl),  discard_l_pl)
         discard_plA%vh(i,1:nsppl)  = pack(symba_plA%helio%swiftest%vh(i,1:npl),  discard_l_pl)
         discard_plA%xb(i,1:nsppl)  = pack(symba_plA%helio%swiftest%xb(i,1:npl),  discard_l_pl)
         discard_plA%vb(i,1:nsppl)  = pack(symba_plA%helio%swiftest%vb(i,1:npl),  discard_l_pl)
         discard_plA%ip(i,1:nsppl)  = pack(symba_plA%helio%swiftest%ip(i,1:npl),  discard_l_pl)
         discard_plA%rot(i,1:nsppl) = pack(symba_plA%helio%swiftest%rot(i,1:npl), discard_l_pl)
      end do

      ! Pack kept bodies down 
      symba_plA%helio%swiftest%name(1:nkpl)   = pack(symba_plA%helio%swiftest%name(1:npl),   .not. discard_l_pl)
      symba_plA%helio%swiftest%status(1:nkpl) = pack(symba_plA%helio%swiftest%status(1:npl), .not. discard_l_pl)
      symba_plA%helio%swiftest%mass(1:nkpl)   = pack(symba_plA%helio%swiftest%mass(1:npl),   .not. discard_l_pl)
      symba_plA%helio%swiftest%radius(1:nkpl) = pack(symba_plA%helio%swiftest%radius(1:npl), .not. discard_l_pl)
      symba_plA%helio%swiftest%rhill(1:nkpl)  = pack(symba_plA%helio%swiftest%rhill(1:npl),  .not. discard_l_pl)
      do i = 1, NDIM
         symba_plA%helio%swiftest%xh(i,1:nkpl)  = pack(symba_plA%helio%swiftest%xh(i,1:npl),  .not. discard_l_pl)
         symba_plA%helio%swiftest%vh(i,1:nkpl)  = pack(symba_plA%helio%swiftest%vh(i,1:npl),  .not. discard_l_pl)
         symba_plA%helio%swiftest%xb(i,1:nkpl)  = pack(symba_plA%helio%swiftest%xb(i,1:npl),  .not. discard_l_pl)
         symba_plA%helio%swiftest%vb(i,1:nkpl)  = pack(symba_plA%helio%swiftest%vb(i,1:npl),  .not. discard_l_pl)
         symba_plA%helio%swiftest%ip(i,1:nkpl)  = pack(symba_plA%helio%swiftest%ip(i,1:npl),  .not. discard_l_pl)
         symba_plA%helio%swiftest%rot(i,1:nkpl) = pack(symba_plA%helio%swiftest%rot(i,1:npl), .not. discard_l_pl)
         symba_plA%helio%ah(i,1:nkpl)           = pack(symba_plA%helio%ah(i,1:npl),           .not. discard_l_pl)
      end do

      if (nkpl + nmergeadd > npl) call util_resize_pl(symba_plA, nkpl+nmergeadd, npl)
      npl = nkpl  + nmergeadd

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

      symba_plA%helio%swiftest%nbody = npl
   end if 

   if (ldiscard_tp) then 
      nktp = 0
      nsptp = 0  

      discard_l_tp(1:ntp) = (symba_tpA%helio%swiftest%status(1:ntp) /= ACTIVE)
      nsptp = count(discard_l_tp)
      nktp = ntp - nsptp

      call discard_tpA%alloc(nsptp) 

      discard_tpA%name(1:nsptp)   = pack(symba_tpA%helio%swiftest%name(1:ntp),   discard_l_tp)
      discard_tpA%status(1:nsptp) = pack(symba_tpA%helio%swiftest%status(1:ntp), discard_l_tp)
      discard_tpA%isperi(1:nsptp) = pack(symba_tpA%helio%swiftest%isperi(1:ntp), discard_l_tp)
      discard_tpA%peri(1:nsptp)   = pack(symba_tpA%helio%swiftest%peri(1:ntp),   discard_l_tp)
      discard_tpA%atp(1:nsptp)    = pack(symba_tpA%helio%swiftest%atp(1:ntp),    discard_l_tp)
      do i = 1, NDIM
         discard_tpA%xh(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%xh(i,1:ntp),   discard_l_tp)
         discard_tpA%vh(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%vh(i,1:ntp),   discard_l_tp)
         discard_tpA%xb(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%xb(i,1:ntp),   discard_l_tp)
         discard_tpA%vb(i,1:nsptp)   = pack(symba_tpA%helio%swiftest%vb(i,1:ntp),   discard_l_tp)
      end do

      symba_tpA%helio%swiftest%name(1:nktp)   = pack(symba_tpA%helio%swiftest%name(1:ntp),   .not. discard_l_tp)
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
   end if 

end subroutine symba_rearray