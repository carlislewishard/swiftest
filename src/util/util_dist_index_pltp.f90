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
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : l            : length of the distance array
!              : ik           : 
!              : jk
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index(npl, swifter_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_index_pltp(nplm, ntp, symba_tpA)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     use module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_pltp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: nplm, ntp
     type(symba_tp), intent(inout) :: symba_tpA

! Internals
     INTEGER(I8B)              :: i,j,ii,jj,np,counter,ii_end,jj_end

! Executable code
     symba_tpA%helio%swiftest%num_pltp_comparisons = (int(nplm, kind = I8B) - 1_I8B) * int(ntp, kind = I8B) ! number of entries in our distance array

     if (allocated(symba_tpA%helio%swiftest%k_pltp)) deallocate(symba_tpA%helio%swiftest%k_pltp)
     allocate(symba_tpA%helio%swiftest%k_pltp(2, symba_tpA%helio%swiftest%num_pltp_comparisons))

! !!$omp parallel do schedule(static) default(none) &
! !!$omp shared(k_pltp, nplm, ntp) &
! !!$omp private(i, j, counter)
!      do i = 2,nplm
!           counter = (i-2) * ntp + 1
!           do j = 1,ntp
!                k_pltp(1,counter) = i
!                k_pltp(2,counter) = j
!                counter = counter + 1
!           enddo
!      enddo
! !!$omp end parallel do

     np = 500

     counter = 1

     do i = 2,nplm,np
        ii_end = min(i+np-1, nplm)
          do j = 1,ntp,np
            jj_end = min(j+np-1, ntp)
               do ii = i, ii_end
                    do jj = j, jj_end
                         symba_tpA%helio%swiftest%k_pltp(1,counter) = ii
                         symba_tpA%helio%swiftest%k_pltp(2,counter) = jj
                         counter = counter + 1
                    enddo
               enddo
          enddo
     enddo

     RETURN

END SUBROUTINE util_dist_index_pltp
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
