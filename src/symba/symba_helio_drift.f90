!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_drift
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Loop through planets and call Danby drift routine
!
!  Input
!    Arguments : irec       : input recursion level
!                npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_helio_drift(irec, npl, symba_pl1P, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_helio_drift.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_helio_drift(irec, npl, symba_plA, dt)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_drift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)       :: irec, npl
     REAL(DP), INTENT(IN)           :: dt
     TYPE(symba_pl), INTENT(INOUT)  :: symba_plA

! Internals
     INTEGER(I4B)              :: i, ndrift
     REAL(DP)                  :: mu
     integer(I4B), dimension(:), allocatable :: iflag
     real(DP), dimension(:,:), allocatable :: xh, vb
     logical(LGT), dimension(:), allocatable :: dodrift 
     integer(I4B) :: ibad

! Executable code
     mu = symba_plA%helio%swiftest%mass(1)
     allocate(dodrift(npl))
     dodrift(2:npl) = (symba_plA%levelg(2:npl) == irec) .AND. (symba_plA%helio%swiftest%status(2:npl) == ACTIVE) 
     ndrift = count(dodrift(2:npl))
     if (ndrift == 0) return
     do i = 2, npl
       if (dodrift(i)) then
         if (symba_plA%helio%swiftest%xh(1, i) /= symba_plA%helio%swiftest%xh(1, i)) then
            write(*,*) 'planet ',i,' is bad!'
         end if
      end if
   end do
     allocate(xh(NDIM, ndrift))
     allocate(vb(NDIM, ndrift))
     allocate(iflag(ndrift))
     iflag = 0
     do i = 1, NDIM
       xh(i,:) = pack(symba_plA%helio%swiftest%xh(i,2:npl), dodrift(2:npl))
       vb(i,:) = pack(symba_plA%helio%swiftest%vb(i,2:npl), dodrift(2:npl))
     end do
     CALL drift_one(mu, xh(:,1:ndrift), vb(:,1:ndrift), dt, iflag(1:ndrift), ndrift)
     if (any(iflag(:) /= 0)) then
         write(*,*) "symba_helio_drift error"
         WRITE(*, *) " STOPPING "
         CALL util_exit(FAILURE)
     end if
     do i = 1, NDIM
         symba_plA%helio%swiftest%xh(i,2:npl) = unpack(xh(i, 1:ndrift), dodrift(2:npl), symba_plA%helio%swiftest%xh(i,2:npl))
         symba_plA%helio%swiftest%vb(i,2:npl) = unpack(vb(i, 1:ndrift), dodrift(2:npl), symba_plA%helio%swiftest%vb(i,2:npl))
     end do
      

     RETURN

END SUBROUTINE symba_helio_drift
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
