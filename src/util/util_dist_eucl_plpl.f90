!**********************************************************************************************************************************
!
!  Unit Name   : util_dist_eucl_plpl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Calculates the Euclidean distance matrix (but in array form)
!
!  Input
!    Arguments : npl : number of planets
!              : invar : variable that we want to make comparisons between
!              : num_comparisons : number of comparisons to make
!              : k_plpl : matrix to convert linear index k into i,j indices
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : outvar : results of comparisons
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index(invar, num_comparisons, k_plpl_outvar)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_eucl_plpl(invar, outvar, symba_plA)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_eucl_plpl
     USE omp_lib
     IMPLICIT NONE

! Arguments
     REAL(DP),DIMENSION(:,:),INTENT(IN) :: invar
     REAL(DP), DIMENSION(:,:),INTENT(INOUT) :: outvar
     type(symba_pl), intent(inout) :: symba_plA

! Internals
     INTEGER(I8B) :: k
     
! Executable code

   !!$omp parallel do schedule(auto) default(private) &
   !!$omp shared (outvar, invar, num_comparisons, k_plpl)
   do k = 1,symba_plA%helio%swiftest%num_plpl_comparisons
      outvar(:,k) = invar(:,symba_plA%helio%swiftest%k_plpl(2,k)) - invar(:,symba_plA%helio%swiftest%k_plpl(1,k))
   end do
   !!$omp end parallel do

     RETURN

END SUBROUTINE util_dist_eucl_plpl
!**********************************************************************************************************************************
!
!  Author(s)   : Jacob R. Elliott 
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
