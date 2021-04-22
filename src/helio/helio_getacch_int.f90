!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch_int
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl        : number of planets
!                helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_getacch_int(npl, helio_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah3.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_getacch_int(npl, helio_plA)

! Modules
     USE swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch_int
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)       :: npl
     TYPE(helio_pl), INTENT(INOUT)  :: helio_plA

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: rji2, irij3, faci, facj
     !REAL(DP), DIMENSION(NDIM) :: dx
     integer(I8B)              :: k
     real(DP)                  :: dx, dy, dz


! Executable code
     !DO i = 2, npl - 1
     !     DO j = i + 1, npl
     !          dx(:) = helio_plA%swiftest%xh(:,j) - helio_plA%swiftest%xh(:,i)
     !          rji2 = DOT_PRODUCT(dx(:), dx(:))
     !          irij3 = 1.0_DP/(rji2*SQRT(rji2))
     !          faci = helio_plA%swiftest%mass(i)*irij3
     !          facj = helio_plA%swiftest%mass(j)*irij3
     !          helio_plA%ahi(:,i) = helio_plA%ahi(:,i) + facj*dx(:)
     !          helio_plA%ahi(:,i) = helio_plA%ahi(:,j) - faci*dx(:)
     !     END DO
     !END DO


     associate(ahi => helio_plA%ahi, xh => helio_plA%swiftest%xh, &
               mass => helio_plA%swiftest%mass, radius => helio_plA%swiftest%radius, &
               num_plpl_comparisons => helio_plA%swiftest%num_plpl_comparisons, &
               k_plpl => helio_plA%swiftest%k_plpl)
         do k = 1, num_plpl_comparisons
            associate(ik => k_plpl(1, k), jk => k_plpl(2, k))
               dx = xh(1, jk) - xh(1, ik)
               dy = xh(2, jk) - xh(2, ik)
               dz = xh(3, jk) - xh(3, ik)
               rji2 = dx**2 + dy**2 + dz**2
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = mass(ik) * irij3
               facj = mass(jk) * irij3
               ahi(1, ik) = ahi(1, ik) + facj * dx
               ahi(2, ik) = ahi(2, ik) + facj * dy
               ahi(3, ik) = ahi(3, ik) + facj * dz
               ahi(1, jk) = ahi(1, jk) - faci * dx
               ahi(2, jk) = ahi(2, jk) - faci * dy
               ahi(3, jk) = ahi(3, jk) - faci * dz
            end associate
         end do
      end associate

     RETURN

END SUBROUTINE helio_getacch_int
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
