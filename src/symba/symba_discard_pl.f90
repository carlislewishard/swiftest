!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their positions or because they are unbound
!
!  Input
!    Arguments : t           : time
!                npl         : number of planets
!                symba_pl1P  : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!                rmin        : minimum allowed heliocentric radius
!                rmax        : maximum allowed heliocentric radius
!                rmaxu       : maximum allowed heliocentric radius for unbound planets
!                qmin        : minimum allowed pericenter distance
!                qmin_coord  : coordinate frame for qmin
!                qmin_alo    : minimum semimajor axis for qmin
!                qmin_ahi    : maximum semimajor axis for qmin
!                ldiscard    : logical flag indicating that at least one body needs to be discarded this step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl         : number of planets
!                symba_pl1P  : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!                eoffset     : energy offset
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_pl(t, npl, symba_pl1P, symba_pld1P, rmin, rmax, rmaxu, qmin, qmin_coord,
!                                      qmin_alo, qmin_ahi)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_massive5.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_pl(t, npl, ntp, symba_plA, symba_tpA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, ldiscard)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: npl, ntp
     REAL(DP), INTENT(IN)        :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)    :: qmin_coord
     TYPE(symba_pl)              :: symba_plA
     TYPE(symba_tp)              :: symba_tpA
     LOGICAL(LGT), INTENT(INOUT) :: ldiscard

! Internals
     REAL(DP)                  :: msys

! Executable code
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP) .OR. ((qmin >= 0.0_DP) .AND. (qmin_coord == "BARY")))      &
          CALL coord_h2b(npl, symba_plA%helio%swiftest, msys)
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP))                                                           &
          CALL symba_discard_sun_pl(t, npl, ntp, msys, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, rmin, rmax, rmaxu, ldiscard)
     IF (qmin >= 0.0_DP) CALL symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscard)
     RETURN 

END SUBROUTINE symba_discard_pl
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
