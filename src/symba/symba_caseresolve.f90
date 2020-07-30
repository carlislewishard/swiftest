!**********************************************************************************************************************************
!
!  Unit Name   : symba_caseresolve
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_caseresolve(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     symba_plA, nplplenc, plplenc_list, regime, plmaxname, tpmaxname, mres, rres, array_index1_child, array_index2_child, &
     m1, m2, rad1, rad2, x1, x2, v1, v2, mtiny)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_caseresolve
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc
     INTEGER(I4B), INTENT(IN)                         :: nplplenc
     INTEGER(I4B), INTENT(INOUT)                      :: plmaxname, tpmaxname, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt, mtiny
     REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
     REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
     REAL(DP), DIMENSION(:), INTENT(IN)            :: vbs
     REAL(DP), DIMENSION(:), INTENT(INOUT)         :: x1, x2, v1, v2
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     INTEGER(I4B), INTENT(IN)                         :: regime
     INTEGER(I4B), DIMENSION(:), INTENT(INOUT)      :: array_index1_child, array_index2_child

! Internals

! Executable code

          SELECT CASE (regime)

          CASE (COLLRESOLVE_REGIME_DISRUPTION)
               CALL symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               CALL symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, symba_plA, nplplenc, &
               plplenc_list, plmaxname, tpmaxname, mres, rres, m1, m2, rad1, &
               rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               CALL symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
               x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_HIT_AND_RUN)
               CALL symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mres, m1, m2, rad1, rad2, x1, x2, v1, v2, mtiny)

          CASE (COLLRESOLVE_REGIME_MERGE)
               CALL symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
               x2, v1, v2)

          CASE DEFAULT 
               WRITE(*,*) "ERROR IN SYMBA_CASERESOLVE, NO REGIME SELECTED"
          END SELECT


RETURN
END SUBROUTINE symba_caseresolve
