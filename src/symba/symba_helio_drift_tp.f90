!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_drift_tp
!  Unit Type   : subroutine
!  Project   : Swiftest
!  Package   : symba
!  Language    : Fortran 90/95
!
!  Description : Loop through test particles and call Danby drift routine
!
!  Input
!    Arguments : irec     : input recursion level
!          ntp      : number of active test particles
!          symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!          mu       : mass of the Sun
!          dt       : time step
!    Terminal  : none
!    File    : none
!
!  Output
!    Arguments : symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : error message
!    File    : none
!
!  Invocation  : CALL symba_helio_drift_tp(irec, ntp, symba_tp1P, mu, dt)
!
!  Notes     : Adapted from Hal Levison's Swift routine symba5_helio_drift.f
!
!**********************************************************************************************************************************
subroutine symba_helio_drift_tp(irec, ntp, symba_tpa, mu, dt)

! Modules
   use swiftest
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => symba_helio_drift_tp
   IMPLICIT NONE

! Arguments
   integer(I4B), intent(in)    :: irec, ntp
   real(DP), intent(in)      :: mu, dt
   type(symba_tp), intent(inout) :: symba_tpa

! Internals
   integer(I4B)          :: i, n
   integer(I4B), dimension(:), allocatable :: iflag
   real(DP), dimension(:, :), allocatable :: xht, vbt
   logical, dimension(:), allocatable :: tpmask

! Executable code
   if (ntp == 0) return
   allocate(tpmask(ntp))
   allocate(xht, source = symba_tpA%helio%swiftest%xh(:,1:ntp))
   allocate(vbt, source = symba_tpA%helio%swiftest%vb(:,1:ntp))
   tpmask(:) = (symba_tpA%levelg(1:ntp) == irec) .and. (symba_tpA%helio%swiftest%status(1:ntp) == ACTIVE)
   n = count(tpmask(:))
   if (n /= ntp) then
      do i = 1, NDIM
         xht(i, :) = pack(xht(i, :), tpmask(:))
         vbt(i, :) = pack(vbt(i, :), tpmask(:))
      end do
   end if
   allocate(iflag(n))
   
   call drift_one(mu, xht(:,1:n), vbt(:,1:n), dt, iflag(:), n)

   if (n /= ntp) then
      do i = 1, NDIM
         xht(i, :) = unpack(xht(i, :), tpmask(:), symba_tpA%helio%swiftest%xh(i, :))
         vbt(i, :) = unpack(vbt(i, :), tpmask(:), symba_tpA%helio%swiftest%vb(i, :))
      end do   
      iflag(:) = unpack(iflag(:), tpmask(:), iflag(:))
   end if 

   do i = 1, ntp
      if (iflag(i) /= 0) then
         symba_tpA%helio%swiftest%status(i) = DISCARDED_DRIFTERR
         write(*, *) "Particle ", symba_tpA%helio%swiftest%name(i), " lost due to error in Danby drift"
      end if
   end do

   return

end subroutine symba_helio_drift_tp
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date      : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State     : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
