!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_interp_eucl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in democratic heliocentric coordinates, calling the recursive
!                subroutine to descend to the appropriate level to handle close encounters
!
!  Input
!    Arguments : 
!                
!                t              : time
!                npl            : number of planets
!                nplm           : number of planets with mass > mtiny
!                ntp            : number of active test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                dt             : time step
!                eoffset        : energy offset (net energy lost in mergers)
!                mtiny          : smallest self-gravitating mass
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
!                eoffset        : energy offset (net energy lost in mergers)
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_step_interp(t, npl, nplm, ntp, symba_pl1P, symba_tp1P, 
!                                       dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd,
!                                       nmergesub, mergeadd_list, mergesub_list )
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_interp.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step_interp_eucl(t, npl, nplm, ntp, symba_plA, symba_tpA,&
   dt, eoffset, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,&
    mergesub_list, param, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_interp_eucl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: npl, nplm, ntp, nplplenc, npltpenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     type(user_input_parameters), intent(in)          :: param
     INTEGER(I8B), INTENT(IN)                         :: num_plpl_comparisons, num_pltp_comparisons
     INTEGER(I4B), DIMENSION(:,:),INTENT(IN) :: k_plpl 
     INTEGER(I4B), DIMENSION(:,:),INTENT(IN) :: k_pltp

! Internals
     INTEGER(I4B)                                 :: i, irec
     REAL(DP)                                     :: dth, msys
     REAL(DP), DIMENSION(NDIM)                    :: ptb, pte
     REAL(DP), DIMENSION(NDIM, npl)               :: xbeg, xend

! Executable code

     dth = 0.5_DP*dt

     CALL coord_vh2vb(npl, symba_plA%helio%swiftest, msys)

     CALL helio_lindrift(npl, symba_plA%helio%swiftest, dth, ptb)
     IF (ntp > 0) THEN
          CALL coord_vh2vb_tp(ntp, symba_tpA%helio%swiftest, -ptb)
          CALL helio_lindrift_tp(ntp, symba_tpA%helio%swiftest, dth, ptb) 
          DO i = 2, npl
               xbeg(:, i) = symba_plA%helio%swiftest%xh(:,i)
          END DO
     END IF

     CALL symba_getacch_eucl(param%lextra_force, t, npl, nplm, symba_plA, param%j2rp2, param%j4rp4, nplplenc, plplenc_list, &
          num_plpl_comparisons, k_plpl)
     IF (ntp > 0) CALL symba_getacch_tp_eucl(param%lextra_force, t, npl, nplm, ntp, symba_plA, symba_tpA, xbeg, param%j2rp2,&
          param%j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)

     CALL helio_kickvb(npl, symba_plA%helio, dth)
     IF (ntp > 0) CALL helio_kickvb_tp(ntp, symba_tpA%helio, dth)
     irec = -1

     CALL symba_helio_drift(irec, npl, symba_plA, dt)
     IF (ntp > 0) CALL symba_helio_drift_tp(irec, ntp, symba_tpA, symba_plA%helio%swiftest%mass(1), dt)
     irec = 0

     CALL symba_step_recur(t, irec, npl, nplm, ntp, symba_plA, symba_tpA, dt, eoffset, nplplenc, npltpenc,              &
          plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, param)

     IF (ntp > 0) THEN
          DO i = 2, npl
               xend(:, i) = symba_plA%helio%swiftest%xh(:,i)
          END DO
     END IF
     CALL symba_getacch_eucl(param%lextra_force, t+dt, npl, nplm, symba_plA, param%j2rp2, param%j4rp4, nplplenc, plplenc_list, &
          num_plpl_comparisons, k_plpl)
     IF (ntp > 0) CALL symba_getacch_tp_eucl(param%lextra_force, t+dt, npl, nplm, ntp, symba_plA, symba_tpA, xend, &
          param%j2rp2,param%j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)
     CALL helio_kickvb(npl, symba_plA%helio, dth)
     IF (ntp > 0) CALL helio_kickvb_tp(ntp, symba_tpA%helio, dth)
     CALL coord_vb2vh(npl, symba_plA%helio%swiftest)
     CALL helio_lindrift(npl, symba_plA%helio%swiftest, dth, pte)
     IF (ntp > 0) THEN
          CALL coord_vb2vh_tp(ntp, symba_tpA%helio%swiftest, -pte)
          CALL helio_lindrift_tp(ntp, symba_tpA%helio%swiftest, dth, pte)
     END IF

     RETURN

END SUBROUTINE symba_step_interp_eucl
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
