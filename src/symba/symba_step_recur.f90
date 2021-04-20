!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_recur
!  Unit Type   : recursive subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current
!                recursion level, if applicable, and descend to the next deeper level if necessary
!
!  Input
!    Arguments : 
!                t              : time
!                ireci          : input recursion level
!                npl            : number of planets
!                nplm           : number of planets with mass > mtiny
!                ntp            : number of active test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                dt0            : time step (primary time step for overall integration)
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : warning message
!    File      : none
!
!  Invocation  : CALL symba_step_recur(t, ireci, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, nplplenc, npltpenc,
!                                      plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_recur.F
!
!**********************************************************************************************************************************
RECURSIVE SUBROUTINE symba_step_recur(t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, &
   nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, param)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_recur
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: ireci, npl, nplm, ntp, nplplenc, npltpenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt0
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     type(user_input_parameters), intent(inout)       :: param

! Internals
     LOGICAL(LGT)              :: lencounter
     INTEGER(I4B)              :: i, j, irecp, icflg, index_i, index_j, index_pl, index_tp
     REAL(DP)                  :: dtl, dth, sgn, mtiny, rji2, rlim2
     REAL(DP), DIMENSION(NDIM) :: xr, vr, vbs

! Executable code

     mtiny = param%mtiny
     dtl = dt0 / (NTENC**ireci)
     dth = 0.5_DP * dtl
     IF (dtl / dt0 < VSMALL) THEN
          WRITE(*, *) "SWIFTEST Warning:"
          WRITE(*, *) "   In symba_step_recur, local time step is too small"
          WRITE(*, *) "   Roundoff error will be important!"
          call util_exit(FAILURE)
     END IF
     irecp = ireci + 1

     !If we are at the topmost recursion level 
     IF (ireci == 0) THEN
          icflg = 0
          !$omp parallel do schedule (auto) default(private) &
          !$omp shared(symba_plA, plplenc_list, nplplenc, irecp, icflg, ireci, dtl)
          DO i = 1, nplplenc
               IF ((plplenc_list%status(i) == ACTIVE) .AND. (plplenc_list%level(i) == ireci)) THEN
                    index_i  = plplenc_list%index1(i)
                    index_j  = plplenc_list%index2(i)
                    xr(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
                    vr(:) = symba_plA%helio%swiftest%vb(:,index_j) - symba_plA%helio%swiftest%vb(:,index_i)
                    CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(index_i),     &  
                         symba_plA%helio%swiftest%rhill(index_j), dtl, irecp, lencounter,                  &
                         plplenc_list%lvdotr(i))
                    !If two bodies are encountering 
                    IF (lencounter) THEN 
                         rlim2 = (symba_plA%helio%swiftest%radius(index_i) + symba_plA%helio%swiftest%radius(index_j))**2
                         rji2 = DOT_PRODUCT(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                         if (rji2 > rlim2) then
                           !$omp critical
                           icflg = 1
                           symba_plA%levelg(index_i) = irecp
                           symba_plA%levelm(index_i) = MAX(irecp, symba_plA%levelm(index_i))
                           symba_plA%levelg(index_j) = irecp
                           symba_plA%levelm(index_j) = MAX(irecp, symba_plA%levelm(index_j))
                           plplenc_list%level(i) = irecp
                           !$omp end critical
                         end if
                    END IF
               END IF
          END DO
          !$omp end parallel do

          !$omp parallel do schedule (auto) default(private) if(npltpenc > omp_get_max_threads())  &
          !$omp shared(symba_plA, symba_tpA, pltpenc_list, npltpenc, irecp, icflg, ireci, dtl)
          DO i = 1, npltpenc
               IF ((pltpenc_list%status(i) == ACTIVE) .AND. (pltpenc_list%level(i) == ireci)) THEN
                    index_pl  = pltpenc_list%indexpl(i)
                    index_tp  = pltpenc_list%indextp(i)
                    
                    xr(:) = symba_tpA%helio%swiftest%xh(:,index_tp) - symba_plA%helio%swiftest%xh(:,index_pl)
                    vr(:) = symba_tpA%helio%swiftest%vb(:,index_tp) - symba_plA%helio%swiftest%vb(:,index_pl)
                    CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(index_pl), 0.0_DP,   &
                         dtl, irecp, lencounter, pltpenc_list%lvdotr(i))
                    IF (lencounter) THEN
                        rlim2 = symba_plA%helio%swiftest%radius(index_pl)**2
                        rji2 = DOT_PRODUCT(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                        if (rji2 > rlim2) then
                           !$omp critical
                           icflg = 1
                           symba_plA%levelg(index_pl) = irecp
                           symba_plA%levelm(index_pl) = MAX(irecp, symba_plA%levelm(index_pl))
                           symba_tpA%levelg(index_tp) = irecp
                           symba_tpA%levelm(index_tp) = MAX(irecp, symba_tpA%levelm(index_tp))
                           pltpenc_list%level(i) = irecp
                           !$omp end critical
                        end if
                    END IF
               END IF
          END DO
          !$omp end parallel do
          lencounter = (icflg == 1)
          sgn = 1.0_DP
          CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA)
          CALL symba_helio_drift(ireci, npl, symba_plA, dtl)
          IF (ntp > 0) CALL symba_helio_drift_tp(ireci, ntp, symba_tpA, symba_plA%helio%swiftest%mass(1), dtl)
          IF (lencounter) CALL symba_step_recur(t, irecp, npl, nplm, ntp, symba_plA, symba_tpA, dt0, &
                      nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, param)
          sgn = 1.0_DP
          CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA) 
          IF (param%lclose) THEN
               vbs(:) = symba_plA%helio%swiftest%vb(:,1)
               DO i = 1, nplplenc
                    index_i  = plplenc_list%index1(i) 
                    index_j  = plplenc_list%index2(i)
                    IF (((plplenc_list%status(i) == ACTIVE) .AND.                                                                 &
                        (symba_plA%levelg(index_i) >= ireci) .AND.                                                              &
                        (symba_plA%levelg(index_j) >= ireci))) THEN 
                        CALL symba_merge_pl (t, dtl, i, nmergesub, mergesub_list, npl, symba_plA, plplenc_list, param)
                     END IF
               END DO
               DO i = 1, npltpenc
                    index_pl  = pltpenc_list%indexpl(i) 
                    index_tp  = pltpenc_list%indextp(i) 
                    IF ((pltpenc_list%status(i) == ACTIVE) .AND.                                          &
                        (symba_plA%levelg(index_pl) >= ireci) .AND.                                       &
                        (symba_tpA%levelg(index_tp) >= ireci)) THEN                                          
                         CALL symba_merge_tp(t, dtl, i, pltpenc_list, vbs, param%encounter_file, symba_plA, symba_tpA)                    !check later 
                    END IF
               END DO
          END IF
          DO i = 1, nplplenc
               index_i  = plplenc_list%index1(i) 
               index_j  = plplenc_list%index2(i) 
               IF (symba_plA%levelg(index_i) == irecp) symba_plA%levelg(index_i) = ireci
               IF (symba_plA%levelg(index_j) == irecp) symba_plA%levelg(index_j) = ireci
               IF (plplenc_list%level(i) == irecp) plplenc_list%level(i) = ireci
          END DO
          DO i = 1, npltpenc
               index_pl  = pltpenc_list%indexpl(i) 
               index_tp  = pltpenc_list%indextp(i) 
               IF (symba_plA%levelg(index_pl) == irecp) symba_plA%levelg(index_pl) = ireci
               IF (symba_tpA%levelg(index_tp) == irecp) symba_tpA%levelg(index_tp) = ireci
               IF (pltpenc_list%level(i) == irecp) pltpenc_list%level(i) = ireci
          END DO
     ELSE
          DO j = 1, NTENC
               icflg = 0
               DO i = 1, nplplenc
                    IF ((plplenc_list%status(i) == ACTIVE) .AND. (plplenc_list%level(i) == ireci)) THEN
                         index_i  = plplenc_list%index1(i) 
                         index_j  = plplenc_list%index2(i) 
                         xr(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
                         vr(:) = symba_plA%helio%swiftest%vb(:,index_j) - symba_plA%helio%swiftest%vb(:,index_i)
                         CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(index_i),    &
                              symba_plA%helio%swiftest%rhill(index_j), dtl, irecp, lencounter,             &
                              plplenc_list%lvdotr(i))
                         IF (lencounter) THEN
                           rlim2 = (symba_plA%helio%swiftest%radius(index_i) + symba_plA%helio%swiftest%radius(index_j))**2
                           rji2 = DOT_PRODUCT(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                           if (rji2 > rlim2) then
                              icflg = 1
                              symba_plA%levelg(index_i) = irecp
                              symba_plA%levelm(index_i) = MAX(irecp, symba_plA%levelm(index_i))
                              symba_plA%levelg(index_j) = irecp
                              symba_plA%levelm(index_j) = MAX(irecp, symba_plA%levelm(index_j))
                              plplenc_list%level(i) = irecp
                           end if
                         END IF
                    END IF
               END DO
               DO i = 1, npltpenc
                    IF ((pltpenc_list%status(i) == ACTIVE) .AND. (pltpenc_list%level(i) == ireci)) THEN
                         index_pl  = pltpenc_list%indexpl(i) 
                         index_tp  = pltpenc_list%indextp(i) 
                         xr(:) = symba_tpA%helio%swiftest%xh(:,index_tp) - symba_plA%helio%swiftest%xh(:,index_pl)
                         vr(:) = symba_tpA%helio%swiftest%vb(:,index_tp)  - symba_plA%helio%swiftest%vb(:,index_pl) 
                         CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(index_pl), 0.0_DP, &     
                              dtl, irecp, lencounter, pltpenc_list%lvdotr(i))
                         IF (lencounter) THEN
                           rlim2 = symba_plA%helio%swiftest%radius(index_pl)**2
                           rji2 = DOT_PRODUCT(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                           if (rji2 > rlim2) then
                              icflg = 1
                              symba_plA%levelg(index_pl) = irecp
                              symba_plA%levelm(index_pl) = MAX(irecp, symba_plA%levelm(index_pl))
                              symba_tpA%levelg(index_tp) = irecp
                              symba_tpA%levelm(index_tp) = MAX(irecp, symba_tpA%levelm(index_tp))
                              pltpenc_list%level(i) = irecp
                           end if
                         END IF
                    END IF
               END DO
               lencounter = (icflg == 1)
               sgn = 1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA) 
               sgn = -1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA)
               CALL symba_helio_drift(ireci, npl, symba_plA, dtl)
               IF (ntp > 0) CALL symba_helio_drift_tp(ireci, ntp, symba_tpA, symba_plA%helio%swiftest%mass(1), dtl)
               IF (lencounter) CALL symba_step_recur(t, irecp, npl, nplm, ntp, symba_plA, symba_tpA, dt0,  &
                    nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, &
                    mergesub_list, param)
               sgn = 1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA) 
               sgn = -1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_plA, symba_tpA)
               IF (param%lclose) THEN
                    vbs(:) = symba_plA%helio%swiftest%vb(:,1)
                    DO i = 1, nplplenc
                         index_i  = plplenc_list%index1(i) 
                         index_j  = plplenc_list%index2(i) 
                         IF ((plplenc_list%status(i) == ACTIVE) .AND.   &
                             (symba_plA%levelg(index_i) >= ireci) .AND. &
                             (symba_plA%levelg(index_j) >= ireci))  THEN    
                             CALL symba_merge_pl(t, dtl, i, nmergesub, mergesub_list, npl, symba_plA, plplenc_list, param)
                         END IF
                    END DO
                    DO i = 1, npltpenc
                         index_pl  = pltpenc_list%indexpl(i) 
                         index_tp  = pltpenc_list%indextp(i) 
                         IF ((pltpenc_list%status(i) == ACTIVE) .AND.    &
                             (symba_plA%levelg(index_pl) >= ireci) .AND. &
                             (symba_tpA%levelg(index_tp) >= ireci))      &
                              CALL symba_merge_tp(t, dtl, i, pltpenc_list, vbs, param%encounter_file, symba_plA, symba_tpA)                !check that later
                    END DO
               END IF
               DO i = 1, nplplenc
                    index_i  = plplenc_list%index1(i) 
                    index_j  = plplenc_list%index2(i) 
                    IF (symba_plA%levelg(index_i) == irecp) symba_plA%levelg(index_i) = ireci
                    IF (symba_plA%levelg(index_j) == irecp) symba_plA%levelg(index_j) = ireci
                    IF (plplenc_list%level(i) == irecp) plplenc_list%level(i) = ireci
               END DO
               DO i = 1, npltpenc
                    index_pl  = pltpenc_list%indexpl(i) 
                    index_tp  = pltpenc_list%indextp(i) 
                    IF (symba_plA%levelg(index_pl) == irecp) symba_plA%levelg(index_pl) = ireci
                    IF (symba_tpA%levelg(index_tp) == irecp) symba_tpA%levelg(index_tp) = ireci
                    IF (pltpenc_list%level(i) == irecp) pltpenc_list%level(i) = ireci
               END DO
          END DO
     END IF

     RETURN

END SUBROUTINE symba_step_recur
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
