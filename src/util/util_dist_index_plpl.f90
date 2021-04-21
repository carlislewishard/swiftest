!**********************************************************************************************************************************
!
!  Unit Name   : util_dist_index_plpl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Turns i,j indices into k index for use in the Euclidean distance matrix
!
!  Input
!    Arguments : npl          : number of planets
!              : nplm         : number of planets above mtiny
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : num_comparisons     : length of the distance array
!              : k_plpl              : matrix of index conversions (linear index to i, j indices)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index_plpl(npl, nplm, num_comparison, k_plpl)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_index_plpl(npl, nplm, num_comparisons, k_plpl)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_plpl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, nplm
     INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: k_plpl
     INTEGER(I8B), INTENT(OUT) :: num_comparisons

! Internals
     INTEGER(I8B)              :: i,j, counter, npl8, nplm8

     npl8 = int(npl, kind = I8B)
     nplm8 = int(nplm, kind = I8B)
! Executable code
     num_comparisons = ((npl8 - 1_I8B) * (npl8 - 2_I8B) / 2_I8B) - & ! number of entries in a strict lower triangle, nplm x npl, minus first column
                       ((npl8 - nplm8 - 1_I8B) * ((npl8 - nplm8 - 1_I8B) + 1_I8B) / 2_I8B)
     if (allocated(k_plpl)) deallocate(k_plpl) 
     allocate(k_plpl(2, num_comparisons))
     ! this is a 'fancier' code, but so far i think it runs slower
     ! so leaving it in, but commenting it out
     ! i think it's because of the 'mod' call, but i haven't profiled it yet
     ! don't forget to uncomment the 'k' declaration up top!
     ! allocate(k(num_comparisons))

     ! m = ceiling(sqrt(2. * num_comparisons))

     ! k = (/(i, i=1,num_comparisons, 1)/)

     ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
     ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2

     ! brute force the index creation

     do i = 2,nplm
          counter = (i - 2) * npl - i * (i - 1) / 2 + 2
          do j = i+1,npl
               k_plpl(1,counter) = i
               k_plpl(2,counter) = j
               counter = counter + 1
          enddo
     enddo

     RETURN

END SUBROUTINE util_dist_index_plpl
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
