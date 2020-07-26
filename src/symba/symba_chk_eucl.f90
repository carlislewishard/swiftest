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
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plpl_encounters
     LOGICAL(LGT), DIMENSION(:), ALLOCATABLE, intent(inout) :: lvdotr
     INTEGER(I4B), INTENT(OUT)                              :: nplplenc

! Internals
     logical, dimension(num_plpl_comparisons) :: lencounter, loc_lvdotr 
     integer(I8B), dimension(num_plpl_comparisons) :: indnum
     ! LOGICAL(LGT) :: iflag lvdotr_flag
     REAL(DP)     :: rcrit, r2crit, vdotr, r2, v2, tmin, r2min, term2, rcritmax, r2critmax
     INTEGER(I8B) :: i, k
     REAL(DP), DIMENSION(NDIM):: xr, vr

    

! Executable code

     nplplenc = 0

     term2 = RHSCALE*RSHELL**0

     rcritmax = (symba_plA%helio%swiftest%rhill(2) + symba_plA%helio%swiftest%rhill(3)) * term2
     r2critmax = rcritmax * rcritmax

      !$omp parallel do default(private) schedule(static) &
      !$omp num_threads(min(omp_get_max_threads(),ceiling(num_plpl_comparisons/10000.))) &
      !$omp shared(num_plpl_comparisons, loc_lvdotr, lencounter, k_plpl, dt, term2, r2critmax, symba_plA, indnum)
      do k = 1, num_plpl_comparisons
         indnum(k) = k
         lencounter(k) = .false.
         xr(:) = symba_plA%helio%swiftest%xh(:,k_plpl(2,k)) - symba_plA%helio%swiftest%xh(:,k_plpl(1,k))

         r2 = dot_product(xr(:), xr(:)) 
         if (r2 < r2critmax) then
            rcrit = (symba_plA%helio%swiftest%rhill(k_plpl(2,k)) + symba_plA%helio%swiftest%rhill(k_plpl(1,k))) * term2
            r2crit = rcrit*rcrit 

            vr(:) = symba_plA%helio%swiftest%vh(:,k_plpl(2,k)) - symba_plA%helio%swiftest%vh(:,k_plpl(1,k))

            vdotr = dot_product(vr(:), xr(:))

            loc_lvdotr(k) = (vdotr < 0.0_DP)

            if (r2 < r2crit) then
               lencounter(k) = .true.
            else
               if (vdotr < 0.0_dp) then
                  v2 = dot_product(vr(:), vr(:))
                  tmin = -vdotr /  v2
                  if (tmin < dt) then
                     r2min = r2 - vdotr * vdotr / v2
                  else
                     r2min = r2 + 2 * vdotr * dt + v2 * dt * dt
                  end if
                  r2min = min(r2min, r2)
                  if (r2min <= r2crit) lencounter(k) = .true.
               end if
            end if
         endif
      end do
     !$omp end parallel do
      nplplenc = count(lencounter(:))
      if (nplplenc > 0) then
         allocate(plpl_encounters(nplplenc))
         allocate(lvdotr(nplplenc))
         plpl_encounters(:) = pack(indnum(:), lencounter(:))
         lvdotr(:) = pack(loc_lvdotr(:), lencounter(:))
      end if
      return

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
