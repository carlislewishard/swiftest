!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_eucl
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
!                lextra_force   : logical flag indicating whether to include user-supplied accelerations
!                lclose         : logical flag indicating whether to check for mergers
!                t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                mtiny          : smallest self-gravitating mass
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
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
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_step(lfirst, lextra_force, lclose, t, npl, ntp, symba_pl1P, symba_tp1P, j2rp2, j4rp4,
!                                dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,
!                                mergesub_list, mtiny, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step_eucl(t,dt,param,npl, ntp,symba_plA, symba_tpA,       &
   nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, &
   mergeadd_list, mergesub_list)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_helio
     USE module_symba
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_eucl
     IMPLICIT NONE

! Arguments
     TYPE(user_input_parameters), INTENT(INOUT)       :: param        ! Derived type containing user defined parameters 
     INTEGER(I4B), INTENT(IN)                         :: npl, ntp
     INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT)              :: lencounter
     INTEGER(I4B)              :: irec, nplm, ipl, ipl1, ipl2
     INTEGER(I8B)              :: i, k, counter
     INTEGER(I8B), DIMENSION(:), ALLOCATABLE :: pltp_encounters_indices
     LOGICAL, SAVE             :: lfirst = .true.
     
     LOGICAL(LGT), ALLOCATABLE, DIMENSION(:) :: pltp_lencounters
     LOGICAL(lgt), ALLOCATABLE, DIMENSION(:) :: plpl_lvdotr, pltp_lvdotr
     
! Executable code
     call symba_step_reset(npl, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list)
     nplplenc = 0
     npltpenc = 0
     irec = 0
     nplm = count(symba_plA%helio%swiftest%mass(1:npl) > param%mtiny)

     CALL symba_chk_eucl(npl, irec, symba_plA, dt, plplenc_list, nplplenc)
     
     if(ntp>0)then
         allocate(pltp_lencounters(symba_tpA%helio%swiftest%num_pltp_comparisons))
         allocate(pltp_lvdotr(symba_tpA%helio%swiftest%num_pltp_comparisons))
         pltp_lencounters = .false.
         pltp_lvdotr = .false.

          CALL symba_chk_eucl_pltp(symba_plA, symba_tpA, dt, pltp_lencounters, pltp_lvdotr, npltpenc)
     
          ! npltpenc = count(pltp_encounters > 0)
          ! print *,'step npltpenc: ',npltpenc
          if(npltpenc>0)then

               allocate(pltp_encounters_indices(npltpenc))

               counter = 1
               do k = 1,symba_tpA%helio%swiftest%num_pltp_comparisons
                    if(pltp_lencounters(k))then
                         pltp_encounters_indices(counter) = k
                         counter = counter + 1
                    endif
               enddo

               symba_plA%ntpenc(symba_tpA%helio%swiftest%k_pltp(1,pltp_encounters_indices(:))) = symba_plA%ntpenc(symba_tpA%helio%swiftest%k_pltp(1,pltp_encounters_indices(:))) + 1
               symba_tpA%nplenc(symba_tpA%helio%swiftest%k_pltp(2,pltp_encounters_indices(:))) = symba_tpA%nplenc(symba_tpA%helio%swiftest%k_pltp(2,pltp_encounters_indices(:))) + 1

               pltpenc_list%status(1:npltpenc) = ACTIVE
               pltpenc_list%lvdotr(1:npltpenc) = pltp_lvdotr(pltp_encounters_indices(:))
               pltpenc_list%level(1:npltpenc) = 0
               pltpenc_list%indexpl(1:npltpenc) = symba_tpA%helio%swiftest%k_pltp(1,pltp_encounters_indices(:))
               pltpenc_list%indextp(1:npltpenc) = symba_tpA%helio%swiftest%k_pltp(1,pltp_encounters_indices(:))

               deallocate(pltp_encounters_indices)
          endif

          deallocate(pltp_lencounters, pltp_lvdotr)
     endif
     
! END OF THINGS THAT NEED TO BE CHANGED IN THE TREE

     ! flag to see if there was an encounter
     lencounter = ((nplplenc > 0) .OR. (npltpenc > 0))

     IF (lencounter) THEN ! if there was an encounter, we need to enter symba_step_interp to see if we need recursion
          CALL symba_step_interp_eucl(t, npl, nplm, ntp, symba_plA, symba_tpA, dt, nplplenc, npltpenc, &
                  plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, param)
          lfirst = .TRUE.
     ELSE ! otherwise we can just advance the particles
         CALL symba_step_helio(lfirst, param%lextra_force, t, npl, nplm, ntp,&
                                 symba_plA%helio, symba_tpA%helio, &
                                 param%j2rp2, param%j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE symba_step_eucl
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
