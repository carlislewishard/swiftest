!**********************************************************************************************************************************
!
!  Unit Name   : util_crossproduct
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Calculates cross product of two arrays
!
!
!  Invocation  : CALL util_crossproduct(symba_plA, npl)
!
!**********************************************************************************************************************************
SUBROUTINE util_crossproduct(ar1, ar2, ans)

! Modules
     USE swiftest_globals
     USE module_interfaces, EXCEPT_THIS_ONE => util_crossproduct
     IMPLICIT NONE

! Arguments
     real(DP),dimension(:),intent(in)  :: ar1,ar2
     real(DP),dimension(:),intent(out) :: ans

! Internals


! Executable code


     ans(1) = ar1(2) * ar2(3) - ar1(3) * ar2(2)
     ans(2) = ar1(3) * ar2(1) - ar1(1) * ar2(3)
     ans(3) = ar1(1) * ar2(2) - ar1(2) * ar2(1)
     
     RETURN

END SUBROUTINE util_crossproduct
!**********************************************************************************************************************************
!
!  Author(s)   : C.Wishard and J.Pouplin
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
