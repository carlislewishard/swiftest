!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_peri_pl
!  Unit Type   : subroutine
!  Project   : Swiftest
!  Package   : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their pericenter distances
!
!  Input
!    Arguments : t      : time
!          npl      : number of planets
!          symba_pl1P : pointer to head of SyMBA planet structure linked-list
!          msys     : total system mass
!          qmin     : minimum allowed pericenter distance
!          qmin_alo   : minimum semimajor axis for qmin
!          qmin_ahi   : maximum semimajor axis for qmin
!          qmin_coord : coordinate frame for qmin
!          ldiscard  : logical flag indicating whether any planets are discarded
!    Terminal  : none
!    File    : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!          ldiscard  : logical flag indicating whether any planets are discarded
!    Terminal  : status message
!    File    : none
!
!  Invocation  : CALL symba_discard_peri_pl(t, npl, symba_pl1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscard)
!
!  Notes     : Adapted from Hal Levison's Swift routine discard_mass_peri.f
!
!**********************************************************************************************************************************
subroutine symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscard)

! modules
   use swiftest
   use module_helio
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_discard_peri_pl
   implicit none

! arguments
   logical(LGT), intent(inout)   :: ldiscard
   integer(I4B), intent(in)    :: npl
   real(dp), intent(in)      :: t, msys, qmin, qmin_alo, qmin_ahi
   character(*), intent(in)    :: qmin_coord
   type(symba_pl), intent(inout) :: symba_plA

! internals
   logical(LGT), save      :: lfirst = .true.
   integer(I4B)          :: i


! executable code
   if (lfirst) then
      call symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
      lfirst = .false.
   else
      call symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
      do i = 2, npl
         if (symba_plA%helio%swiftest%status(i) == ACTIVE) then
            if ((symba_plA%isperi(i) == 0) .and. (symba_plA%nplenc(i)== 0)) then
               if ((symba_plA%atp(i) >= qmin_alo) .and. (symba_plA%atp(i) <= qmin_ahi) &
                .and. (symba_plA%peri(i) <= qmin)) then
                  ldiscard = .true.
                  symba_plA%helio%swiftest%status(i) = DISCARDED_PERI
                  write(*, *) "Particle ", symba_plA%helio%swiftest%name(i), &
                   " perihelion distance too small at t = ", t
               end if
            end if
         end if
      end do
   end if

   return

end subroutine symba_discard_peri_pl