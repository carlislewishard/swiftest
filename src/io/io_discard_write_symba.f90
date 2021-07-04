!**********************************************************************************************************************************
!
!  Unit Name   : io_discard_write_symba
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write out information about discarded and merged planets and test particles in SyMBA
!
!  Input
!    Arguments : t             : time
!                mtiny         : smallest self-gravitating mass
!                npl           : number of planets
!                nsppl         : number of spilled planets
!                nsptp         : number of spilled test particles
!                nmergeadd     : number of merged planets to add
!                nmergesub     : number of merged planets to subtract
!                symba_pl1P    : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P   : pointer to head of discard SyMBA planet structure linked-list
!                symba_tpd1P   : pointer to head of discard SyMBA test particle structure linked-list
!                mergeadd_list : array of structures of merged planets to add
!                mergesub_list : array of structures of merged planets to subtract
!                fname         : name of file to write
!                lbig_discard  : logical flag indicating whether to dump planet data with discards
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : discard data to discard file
!
!  Invocation  : CALL io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergeadd, nmergesub, symba_pl1P, symba_pld1P,
!                                            symba_tpd1P, mergeadd_list, mergesub_list, fname, lbig_discard)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_discard_mass.f and io_discard_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergesub, symba_plA, & 
     discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)

! Modules
     USE swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => io_discard_write_symba
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                       :: lbig_discard
     INTEGER(I4B), INTENT(IN)                       :: npl, nsppl, nsptp, nmergesub
     REAL(DP), INTENT(IN)                           :: t, mtiny
     CHARACTER(*), INTENT(IN)                       :: fname
     TYPE(symba_pl), INTENT(INOUT)                  :: symba_plA
     TYPE(swiftest_tp), INTENT(INOUT)               :: discard_tpA
     TYPE(swiftest_pl), INTENT(INOUT)               :: discard_plA
     TYPE(symba_merger), INTENT(INOUT)              :: mergeadd_list, mergesub_list

! Internals
     INTEGER(I4B), PARAMETER   :: LUN = 40
     INTEGER(I4B)              :: i, index, ierr, nplm, nadded

! Executable code
     CALL io_open(LUN, fname, "APPEND", "FORMATTED", ierr)
     IF (ierr /= 0) THEN
          CALL io_open(LUN, fname, "NEW", "FORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTEST Error:"
               WRITE(*, *) "   Unable to open discard output file, ", fname
               CALL util_exit(FAILURE)
          END IF
     END IF
     WRITE(LUN, 100) t, nsppl + nsptp, lbig_discard 
 100 FORMAT(E23.16, 1X, I8, 1X, L1)
     index = 0
     DO i = 1, (nmergesub / 2)
          WRITE(LUN, 200) SUB, mergesub_list%id(i), mergesub_list%status(i)
 200      FORMAT(A, 2(1X, I8))
          WRITE(LUN, 300) mergesub_list%xb(:,i) - symba_plA%helio%swiftest%xb(:,1)
 300      FORMAT(3(E23.16, 1X))
          WRITE(LUN, 300) mergesub_list%vb(:,i) - symba_plA%helio%swiftest%vb(:,1)
          nadded = mergesub_list%nadded(i)

          WRITE(LUN, 200) SUB, mergesub_list%id(i + 1), mergesub_list%status(i + 1)
          WRITE(LUN, 300) mergesub_list%xb(:,i + 1) - symba_plA%helio%swiftest%xb(:,1)
          WRITE(LUN, 300) mergesub_list%vb(:,i + 1) - symba_plA%helio%swiftest%vb(:,1)
          WRITE(LUN, 500) mergesub_list%mass(i), mergesub_list%radius(i)

          DO index = 1, nadded
               WRITE(LUN, 200) ADD, mergeadd_list%id(index), mergeadd_list%status(index)
               WRITE(LUN, 300) mergeadd_list%xb(:,index) - symba_plA%helio%swiftest%xb(:,1)
               WRITE(LUN, 300) mergeadd_list%vb(:,index) - symba_plA%helio%swiftest%vb(:,1)
               WRITE(LUN, 500) mergeadd_list%mass(index), mergeadd_list%radius(index)
          END DO
     END DO

     DO i = 1, (nsppl / 2)
          IF ((discard_plA%status(i) /= MERGED) .AND. (discard_plA%status(i) /= HIT_AND_RUN) .AND. &
             (discard_plA%status(i) /= DISRUPTION) .AND. (discard_plA%status(i) /= SUPERCATASTROPHIC)) THEN
               WRITE(LUN, 200) SUB, discard_plA%id(i), discard_plA%status(i)
               WRITE(LUN, 300) discard_plA%xh(1,i),discard_plA%xh(2,i), discard_plA%xh(3,i)
               WRITE(LUN, 300) discard_plA%vh(1,i),discard_plA%vh(2,i), discard_plA%vh(3,i)
               WRITE(LUN, 500) discard_plA%mass(i), discard_plA%radius(i)
          END IF
     END DO
     DO i = 1, nsptp
          WRITE(LUN, 200) SUB, discard_tpA%id(i), discard_tpA%status(i)
          WRITE(LUN, 300) discard_tpA%xh(1,i),discard_tpA%xh(2,i), discard_tpA%xh(3,i)
          WRITE(LUN, 300) discard_tpA%vh(1,i),discard_tpA%vh(2,i), discard_tpA%vh(3,i)
     END DO
     IF (lbig_discard) THEN
          nplm = 0
          DO i = 1, npl
               IF (symba_plA%helio%swiftest%mass(i) < mtiny) EXIT
               nplm = nplm + 1
          END DO
          IF (nplm > 1) THEN
               WRITE(LUN, 400) nplm
 400           FORMAT(I8)
               DO i = 2, nplm
                    WRITE(LUN, 600) symba_plA%helio%swiftest%id(i), symba_plA%helio%swiftest%mass(i),& 
                     symba_plA%helio%swiftest%radius(i)
 500                FORMAT(2(1X, E23.16))
 600                FORMAT(I8, 2(1X, E23.16))
                    WRITE(LUN, 300) symba_plA%helio%swiftest%xh(:,i)
                    WRITE(LUN, 300) symba_plA%helio%swiftest%vh(:,i)
               END DO
          END IF
     END IF
     CLOSE(LUN)

     RETURN

END SUBROUTINE io_discard_write_symba
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