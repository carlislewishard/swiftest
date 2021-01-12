!**********************************************************************************************************************************
!
!  Unit Name   : symba_step
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in democratic heliocentric coordinates, descending the recursive
!                branch if necessary to handle possible close encounters
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                param%j2rp2          : J2 * R**2 for the Sun
!                param%j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_step(lfirst, param%lextra_force, param%lclose, t, npl, ntp, symba_pl1P, symba_tp1P, param%j2rp2, param%j4rp4,
!                                dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,
!                                mergesub_list, eoffset)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step(t, dt, param, npl, ntp,symba_plA, symba_tpA,       &
               nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub,&
               mergeadd_list, mergesub_list, eoffset)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step
     IMPLICIT NONE

! Arguments
     TYPE(user_input_parameters), INTENT(INOUT)       :: param        ! Derived type containing user defined parameters 
     INTEGER(I4B), INTENT(IN)                         :: npl, ntp
     INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
! Internals
     LOGICAL(LGT)              :: lencounter, lvdotr
     INTEGER(I4B)              :: i, j, irec, nplm
     REAL(DP), DIMENSION(NDIM) :: xr, vr
     LOGICAL, SAVE             :: lfirst = .true.

! Executable code

   symba_plA%index_parent(1:npl) = (/ (i, i=1, npl) /)
   symba_plA%nplenc(:) = 0
   symba_plA%ntpenc(:) = 0
   symba_plA%levelg(:) = -1
   symba_plA%levelm(:) = -1
   symba_plA%index_child(:, :) = 0
   symba_tpA%nplenc(:) = 0 
   symba_tpA%levelg(:) = -1
   symba_tpA%levelm(:) = -1

   nplplenc = 0
   npltpenc = 0
   irec = 0

   nplm = count(symba_plA%helio%swiftest%mass(1:npl) >= param%mtiny)
   if (.not. allocated(plplenc_list%status)) call symba_plplenc_allocate(plplenc_list, 1)
   if (.not. allocated(pltpenc_list%status)) call symba_pltpenc_allocate(pltpenc_list, 1)
   if (.not. allocated(mergeadd_list%status)) call symba_merger_allocate(mergeadd_list, 1)
   if (.not. allocated(mergesub_list%status)) call symba_merger_allocate(mergesub_list, 1)

     DO i = 2, nplm
         !$omp parallel do schedule(static) default(private) &
         !$omp shared(i, npl, nplm, symba_plA, param, dt, irec, plplenc_list, nplplenc)
         DO j = i + 1, npl
               xr(:) = symba_plA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
               vr(:) = symba_plA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
               CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), &
                    symba_plA%helio%swiftest%rhill(j), dt, irec, lencounter, lvdotr)
               IF (lencounter) THEN
                  !$omp critical
                  nplplenc = nplplenc + 1
                  call symba_plplenc_size_check(plplenc_list, nplplenc)
                  plplenc_list%status(nplplenc) = ACTIVE
                  plplenc_list%lvdotr(nplplenc) = lvdotr
                  plplenc_list%level(nplplenc) = irec
                  plplenc_list%index1(nplplenc) = i
                  plplenc_list%index2(nplplenc) = j
                  symba_plA%lmerged(i) = .FALSE.
                  symba_plA%nplenc(i) = symba_plA%nplenc(i) + 1
                  symba_plA%levelg(i) = irec
                  symba_plA%levelm(i) = irec
                  symba_plA%nchild(i) = 0 
                  symba_plA%lmerged(j) = .FALSE.
                  symba_plA%nplenc(j) = symba_plA%nplenc(j) + 1
                  symba_plA%levelg(j) = irec
                  symba_plA%levelm(j) = irec
                  symba_plA%nchild(j) = 0
                  !$omp end critical
               END IF
          END DO
         !$omp end parallel do
          DO j = 1, ntp
               xr(:) = symba_tpA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
               vr(:) = symba_tpA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
               CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), 0.0_DP, dt, irec, lencounter, lvdotr)
               IF (lencounter) THEN
                  npltpenc = npltpenc + 1
                  call symba_pltpenc_size_check(pltpenc_list, npltpenc)
                  symba_plA%ntpenc(i) = symba_plA%ntpenc(i) + 1
                  symba_plA%levelg(i) = irec
                  symba_plA%levelm(i) = irec
                  symba_tpA%nplenc(j) = symba_tpA%nplenc(j) + 1
                  symba_tpA%levelg(j) = irec
                  symba_tpA%levelm(j) = irec
                  pltpenc_list%status(npltpenc) = ACTIVE
                  pltpenc_list%lvdotr(npltpenc) = lvdotr
                  pltpenc_list%level(npltpenc) = irec
                  pltpenc_list%indexpl(npltpenc) = i
                  pltpenc_list%indextp(npltpenc) = j
               END IF
          END DO
     END DO

     lencounter = ((nplplenc > 0) .OR. (npltpenc > 0))
     IF (lencounter) THEN
          CALL symba_step_interp(t, npl, nplm, ntp, symba_plA, symba_tpA, &
               dt, eoffset, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
               nmergesub, mergeadd_list, mergesub_list,  param)
          lfirst = .TRUE.
     ELSE 
          CALL symba_step_helio(lfirst, param%lextra_force, t, npl, nplm, ntp,&
               symba_plA%helio, symba_tpA%helio, &
               param%j2rp2, param%j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE symba_step
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
