!**********************************************************************************************************************************
!
!  Unit Name   : symba_chk_eucl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check for an encounter
!
!  Input
!    Arguments : xr         : position of body 2 relative to body 1
!                vr         : velocity of body 2 relative to body 1
!                rhill1     : Hill sphere radius of body 1
!                rhill2     : Hill sphere radius of body 2
!                dt         : time step
!                irec       : recursion level
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lencounter : logical flag indicating whether there is an encounter this time step
!                lvdotr     : logical flag indicating whether the two bodies are approaching
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_chk.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_chk_eucl(num_plpl_comparisons, k_plpl, symba_plA, dt, plpl_encounters, lvdotr, nplplenc)

! Modules
     USE swiftest
     USE swiftest_globals
     USE swiftest_data_structures
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_chk_eucl
     !$ USE omp_lib
     IMPLICIT NONE

! Arguments
     INTEGER(I8B), INTENT(IN)                               :: num_plpl_comparisons
     INTEGER(I4B), DIMENSION(:,:), INTENT(IN)               :: k_plpl
     TYPE(symba_pl), INTENT(IN)                             :: symba_plA
     REAL(DP), INTENT(IN)                                   :: dt
     INTEGER(I8B), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plpl_encounters
     LOGICAL(LGT), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: lvdotr
     INTEGER(I4B), INTENT(OUT)                              :: nplplenc

! Internals
     ! LOGICAL(LGT) :: iflag lvdotr_flag
     REAL(DP)     :: rcrit, r2crit, vdotr, r2, v2, tmin, r2min, term2, rcritmax, r2critmax
     INTEGER(I8B) :: k, counter
     REAL(DP), DIMENSION(NDIM):: xr, vr
     logical, dimension(num_plpl_comparisons) :: loc_lvdotr, lencounter

! Executable code

     nplplenc = 0

     term2 = RHSCALE*RSHELL**0
     lencounter(:) = .false.
     loc_lvdotr(:) = .false.

     rcritmax = (symba_plA%helio%swiftest%rhill(2) + symba_plA%helio%swiftest%rhill(3)) * term2
     r2critmax = rcritmax * rcritmax

     !$omp parallel do default(private) schedule(static) &
     !$omp num_threads(min(omp_get_max_threads(),ceiling(num_plpl_comparisons/10000.))) &
     !$omp shared(num_plpl_comparisons, loc_lvdotr, lencounter, k_plpl, dt, term2, r2critmax, symba_plA) &
     !$omp reduction(+:nplplenc)
     do k = 1,num_plpl_comparisons
          xr(:) = symba_plA%helio%swiftest%xh(:,k_plpl(2,k)) - symba_plA%helio%swiftest%xh(:,k_plpl(1,k))

          r2 = DOT_PRODUCT(xr(:), xr(:)) 
          if (r2<r2critmax) then
               rcrit = (symba_plA%helio%swiftest%rhill(k_plpl(2,k)) + symba_plA%helio%swiftest%rhill(k_plpl(1,k))) * term2
               r2crit = rcrit*rcrit 

               vr(:) = symba_plA%helio%swiftest%vh(:,k_plpl(2,k)) - symba_plA%helio%swiftest%vh(:,k_plpl(1,k))

               vdotr = DOT_PRODUCT(vr(:), xr(:))

               loc_lvdotr(k) = (vdotr < 0.0_DP) 

               IF (r2 < r2crit) THEN
                    lencounter(k) = .true.
                    nplplenc = nplplenc + 1
               ELSE
                    IF (vdotr < 0.0_DP) THEN
                         v2 = DOT_PRODUCT(vr(:), vr(:))
                         tmin = -vdotr/v2
                         IF (tmin < dt) THEN
                              r2min = r2 - vdotr*vdotr/v2
                         ELSE
                              r2min = r2 + 2.0_DP*vdotr*dt + v2*dt*dt
                         END IF
                         r2min = MIN(r2min, r2)
                         IF (r2min <= r2crit) then
                              lencounter(k) = .true.
                              nplplenc = nplplenc + 1
                         endif
                    END IF
               END IF
          endif
     enddo
   !$omp end parallel do
     if (nplplenc == 0) return
     if (allocated(plpl_encounters)) deallocate(plpl_encounters)
     if (allocated(lvdotr)) deallocate(lvdotr)
     allocate(plpl_encounters(nplplenc))
     allocate(lvdotr(nplplenc))
     counter = 1
     do k = 1,num_plpl_comparisons
          if(lencounter(k))then
               plpl_encounters(counter) = k
               lvdotr(counter) = loc_lvdotr(k)
               counter = counter + 1_I8B
          endif
     enddo
     
     RETURN

END SUBROUTINE symba_chk_eucl
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
