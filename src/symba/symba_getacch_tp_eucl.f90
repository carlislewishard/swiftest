!**********************************************************************************************************************************
!
!  Unit Name   : symba_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of test particles
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                ntp          : number of active test particles
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!                xh           : heliocentric positions of planets at time t
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                npltpenc     : number of planet-test particle encounters
!                pltpenc_list : array of planet-test particle encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_getacch_tp(lextra_force, t, npl, ntp, symba_pl1P, symba_tp1P, xh, j2rp2, j4rp4,
!                                      npltpenc, pltpenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_getacch.f
!
!                Accelerations in an encounter are not included here
!
!**********************************************************************************************************************************
SUBROUTINE symba_getacch_tp_eucl(lextra_force, t, npl, ntp, symba_plA, symba_tpA, xh, j2rp2, j4rp4, npltpenc, pltpenc_list)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch_tp_eucl
     USE omp_lib
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, ntp, npltpenc
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     REAL(DP), DIMENSION(:, :), INTENT(IN)         :: xh
     TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
     TYPE(symba_pltpenc), INTENT(IN)               :: pltpenc_list


! Internals
     INTEGER(I8B)                                 :: i, j, k, index_pl, index_tp
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(:), ALLOCATABLE          :: irh, irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE       :: aobl, xht, aoblt
     REAL(DP), DIMENSION(NDIM, ntp)               :: ah

! Executable code

     ah(:,1:ntp) = 0.0_DP

     ! CALL util_dist_eucl_pltp(symba_plA%helio%swiftest%xh, symba_tpA%helio%swiftest%xh, &
     !      num_pltp_comparisons, k_pltp, dist_pltp_array)

!!$omp parallel do default(shared) schedule(static) &
!!$omp shared(num_pltp_comparisons, symba_plA, symba_tpA, k_pltp) &
!!$omp reduction(-:ah)
     do k = 1,symba_tpA%helio%swiftest%num_pltp_comparisons
          j = symba_tpA%helio%swiftest%k_pltp(2,k)
          IF (symba_tpA%helio%swiftest%status(j) == ACTIVE) THEN
               i = symba_tpA%helio%swiftest%k_pltp(1,k)
               dx(:) = symba_tpA%helio%swiftest%xh(:,symba_tpA%helio%swiftest%k_pltp(2,k)) - symba_plA%helio%swiftest%xh(:,symba_tpA%helio%swiftest%k_pltp(1,k))
               r2 = DOT_PRODUCT(dx(:), dx(:))
               fac = symba_PlA%helio%swiftest%mass(i)/(r2*SQRT(r2))
               ah(:,j) = ah(:,j) - fac*dx(:)
          endif
     enddo
!!$omp end parallel do

     symba_tpA%helio%ah(:,1:ntp) = ah(:,1:ntp)

     !!$omp PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
     !!$omp PRIVATE(i,swifter_plP,helio_tpP,dx,r2,fac) &
     !!$omp SHARED(pltpenc_list,npltpenc)
     DO i = 1, npltpenc
          index_pl = pltpenc_list%indexpl(i)
          index_tp = pltpenc_list%indextp(i)
          IF (symba_tpA%helio%swiftest%status(index_tp) == ACTIVE) THEN
               dx(:) = symba_tpA%helio%swiftest%xh(:,index_tp) - symba_plA%helio%swiftest%xh(:,index_pl)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               fac = symba_plA%helio%swiftest%mass(index_pl)/(r2*SQRT(r2))
               symba_tpA%helio%ah(:,index_tp) = symba_tpA%helio%ah(:,index_tp) + fac*dx(:)
          END IF
     END DO
     !!$omp END PARALLEL DO

     IF (j2rp2 /= 0.0_DP) THEN
         if (npl > 0) then
            ALLOCATE(aobl(NDIM, npl), irh(npl))
         end if
         if (ntp > 0) then
            ALLOCATE(xht(NDIM, ntp), aoblt(NDIM, ntp), irht(ntp))
         end if      
         DO i = 2, npl
           r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
           irh(i) = 1.0_DP/SQRT(r2)
         END DO
         CALL obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, aobl)
         mu = symba_plA%helio%swiftest%mass(1)
         DO i = 1, ntp
            xht(:, i) = symba_tpA%helio%swiftest%xh(:,i) 
            r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
            irht(i) = 1.0_DP/SQRT(r2)
         END DO
         CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
         DO i = 1, ntp
            IF (symba_tpA%helio%swiftest%status(i) == ACTIVE) &
            symba_tpA%helio%ah(:,i) = symba_tpA%helio%ah(:,i) + aoblt(:, i) - aobl(:, 1)
         END DO
         deallocate(aobl, irh, xht, aoblt, irht)
     END IF
     IF (lextra_force) CALL symba_user_getacch_tp(t, ntp, symba_tpA)

     RETURN

END SUBROUTINE symba_getacch_tp_eucl
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
