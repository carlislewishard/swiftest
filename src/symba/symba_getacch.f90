!**********************************************************************************************************************************
!
!  Unit Name   : symba_getacch
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplm         : number of planets with mass > mtiny
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-test particle encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_getacch(lextra_force, t, npl, nplm, symba_pl1P, j2rp2, j4rp4, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_getacch.f
!
!                Accelerations in an encounter are not included here
!
!**********************************************************************************************************************************
SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplplenc
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
     TYPE(symba_plplenc), INTENT(INOUT)            :: plplenc_list

! Internals
     INTEGER(I4B)                                 :: i, j, index_i, index_j
     REAL(DP)                                     :: rji2, irij3, faci, facj, r2, rlim2
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(npl)                     :: irh
     REAL(DP), DIMENSION(NDIM, npl)               :: aobl
     REAL(DP), DIMENSION(NDIM, npl)               :: ahp, ahm


! Executable code
     
     ahp(:,:) = 0.0_DP
     ahm(:,:) = 0.0_DP
     symba_plA%helio%ah(:,:) = 0.0_DP

     !$omp parallel do schedule(static) default(private) &
     !$omp shared(nplm, npl, symba_plA) &
     !$omp reduction(+:ahp) &
     !$omp reduction(-:ahm)
     DO i = 2, nplm
          DO j = i + 1, npl
               IF ((.NOT. symba_plA%lmerged(i)) .OR. (.NOT. symba_plA%lmerged(j)) .OR. &
                   (symba_plA%index_parent(i) /= symba_plA%index_parent(j))) THEN
                    dx(:) = symba_plA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
                    rlim2 = (symba_plA%helio%swiftest%radius(i) + symba_plA%helio%swiftest%radius(j))**2
                    rji2 = DOT_PRODUCT(dx(:), dx(:))
                    if (rji2 > rlim2) then !If false, we likely have recent fragments with coincident positions. 
                                           !  so ignore in this step and let the merge code deal with it in the next
                       irij3 = 1.0_DP / (rji2 * SQRT(rji2))
                       faci = symba_plA%helio%swiftest%mass(i)*irij3
                       facj = symba_plA%helio%swiftest%mass(j)*irij3
                       ahp(:, i) = ahp(:, i) + facj * dx(:)
                       ahm(:, j) = ahm(:, j) - faci * dx(:)
                    end if
               END IF
          END DO
     END DO
     !$omp end parallel do
     symba_plA%helio%ah(:,1:npl) = ahp(:, :) + ahm(:,:)

     DO i = 1, nplplenc
          index_i = plplenc_list%index1(i)
          index_j = plplenc_list%index2(i)
          IF ((.NOT. symba_plA%lmerged(index_i)) .OR. (.NOT. symba_plA%lmerged(index_j))  &
                .OR. (symba_plA%index_parent(index_i) /= symba_plA%index_parent(index_j))) THEN !need to update parent/children
               dx(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
               rlim2 = (symba_plA%helio%swiftest%radius(index_i) + symba_plA%helio%swiftest%radius(index_j))**2
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               if (rji2 > rlim2) then !If false, we likely have recent fragments with coincident positions. 
                                      !  so ignore in this step and let the merge code deal with it in the next
                  irij3 = 1.0_DP / (rji2 * SQRT(rji2))
                  faci = symba_plA%helio%swiftest%mass(index_i) * irij3
                  facj = symba_plA%helio%swiftest%mass(index_j) * irij3
                  symba_plA%helio%ah(:, index_i) = symba_plA%helio%ah(:, index_i) - facj * dx(:)
                  symba_plA%helio%ah(:, index_j) = symba_plA%helio%ah(:, index_j) + faci * dx(:)
               end if
          END IF
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          DO i = 2, npl
               r2 = DOT_PRODUCT(symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%xh(:,i))
               irh(i) = 1.0_DP / SQRT(r2)
          END DO
          CALL obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, aobl)
          DO i = 2, npl
               symba_plA%helio%ah(:,i) = symba_plA%helio%ah(:,i) + aobl(:, i) - aobl(:, 1)
          END DO
     END IF
     IF (lextra_force) CALL symba_user_getacch(t, npl, symba_plA)
     RETURN

END SUBROUTINE symba_getacch
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
