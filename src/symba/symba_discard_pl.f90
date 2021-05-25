subroutine symba_discard_pl(t, npl, ntp, symba_plA, symba_tpA, rmin, rmax, rmaxu, &
                            qmin, qmin_coord, qmin_alo, qmin_ahi, discard_l_pl, discard_stat_list)
   !! author: David A. Minton
   !!
   !! Check to see if planets should be discarded based on their positions or because they are unbound
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_discard_pl.f90
   !! Adapted from Martin Duncan and Hal Levison's Swift routine discard_massive5.f
   use swiftest
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => symba_discard_pl
   implicit none

   ! Arguments
   integer(I4B), intent(inout) :: npl, ntp
   real(DP), intent(in)      :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
   character(*), intent(in)    :: qmin_coord
   type(symba_pl), intent(inout)  :: symba_plA
   type(symba_tp), intent(inout)  :: symba_tpa
   logical, dimension(:), intent(out) :: discard_l_pl
   integer(I4B), dimension(:), allocatable, intent(out) :: discard_stat_list

   ! Internals
   real(DP)            :: msys
   integer(I4B) :: i, ndiscard
   logical :: ldiscard

   ! Executable code
   if ((rmin >= 0.0_DP) .or. (rmax >= 0.0_DP) .or. (rmaxu >= 0.0_DP) .or. ((qmin >= 0.0_DP) .and. (qmin_coord == "bary")))    &
      call coord_h2b(npl, symba_plA%helio%swiftest, msys)
   if ((rmin >= 0.0_DP) .or. (rmax >= 0.0_DP) .or. (rmaxu >= 0.0_DP))                                     &
      call symba_discard_sun_pl(t, npl, ntp, msys, symba_plA%helio%swiftest, symba_tpa%helio%swiftest, rmin, rmax, rmaxu, ldiscard)
   !if (qmin >= 0.0_DP) call symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscard)
   if (.not.ldiscard) return

   ! We need to keep track of the bodies that are discarded via collision or escape from the central body separately from those that collide with each other.
   associate(status => symba_plA%helio%swiftest%status)
      discard_l_pl(1:npl) = (status(1:npl) == DISCARDED_RMIN) .or. (status(1:npl) == DISCARDED_PERI) .or. (status(1:npl) == DISCARDED_RMAX) .or. (status(1:npl) == DISCARDED_RMAXU)
      ndiscard = count(discard_l_pl(1:npl))
      if (ndiscard > 0) then
         if (allocated(discard_stat_list)) deallocate(discard_stat_list)
         allocate(discard_stat_list(ndiscard))
         discard_stat_list = pack(status(1:npl), discard_l_pl)
      end if
   end associate
       
   return 

end subroutine symba_discard_pl