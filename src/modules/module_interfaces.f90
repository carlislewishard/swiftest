!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swiftest
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in swiftest package
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
MODULE module_interfaces
   use swiftest_globals
   use swiftest_data_structures
   use user

     IMPLICIT NONE


     INTERFACE
          SUBROUTINE coord_b2h(npl, swiftest_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_b2h
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_b2h_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b(npl, swiftest_plA, msys)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(OUT)            :: msys
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b_tp(ntp, swiftest_tpA, swiftest_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh(npl, swiftest_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vb2vh
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh_tp(ntp, swiftest_tpA, vs)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(:), INTENT(IN) :: vs
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE coord_vb2vh_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb(npl, swiftest_plA, msys)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: npl
               REAL(DP), INTENT(OUT)           :: msys
               TYPE(swiftest_pl),INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vh2vb
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb_tp(ntp, swiftest_tpA, vs)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(:), INTENT(IN) :: vs
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE coord_vh2vb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE discard(t, dt, npl, ntp, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, qmin,  &
               qmin_alo, qmin_ahi, qmin_coord, lclose, lrhill_present)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_peri(t, npl, ntp, swiftest_plA, swiftest_tpA, msys, qmin, qmin_alo, & 
               qmin_ahi, qmin_coord, lrhill_present)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard_peri
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)             :: iflag
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(:), INTENT(IN) :: dx, dv
               REAL(DP), INTENT(OUT)                 :: r2min
          END SUBROUTINE discard_pl_close
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_sun(t, ntp, msys, swifter_tpA, rmin, rmax, rmaxu)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: t, msys, rmin, rmax, rmaxu
               TYPE(swiftest_tp), INTENT(INOUT) :: swifter_tpA
          END SUBROUTINE discard_sun
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_dan(mu, x0, v0, dt0, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)                :: iflag
               REAL(DP), INTENT(IN)                     :: mu, dt0
               REAL(DP), DIMENSION(:), INTENT(INOUT)    :: x0, v0
          END SUBROUTINE drift_dan
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepmd(dm, es, ec, x, s, c)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dm, es, ec
               REAL(DP), INTENT(OUT) :: x, s, c
          END SUBROUTINE drift_kepmd
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u, s
               REAL(DP), INTENT(OUT) :: f
          END SUBROUTINE drift_kepu_fchk
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_guess(dt, r0, mu, alpha, u, s)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT) :: s
          END SUBROUTINE drift_kepu_guess
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_lag
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_new
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: s
          END SUBROUTINE drift_kepu_p3solve
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_stumpff(x, c0, c1, c2, c3)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(INOUT) :: x
               REAL(DP), INTENT(OUT)   :: c0, c1, c2, c3
          END SUBROUTINE drift_kepu_stumpff
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_one(mu, x, v, dt, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)  :: iflag
               REAL(DP), INTENT(IN)                  :: mu, dt
               REAL(DP), DIMENSION(:), INTENT(INOUT) :: x, v
          END SUBROUTINE drift_one
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift(npl, swiftest_plA, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(IN)             :: dt
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift_tp(ntp, swiftest_tpA, mu, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, dt
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch(lflag, lextra_force, t, npl, helio_plA, j2rp2, j4rp4)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)       :: npl
               REAL(DP), INTENT(IN)           :: t, j2rp2, j4rp4
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int(npl, helio_plA)
               USE swiftest_globals
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: npl
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int_tp(npl, ntp, swiftest_plA, helio_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               TYPE(swiftest_pl), INTENT(INOUT)           :: swiftest_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_getacch_int_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, ntp, helio_plA, helio_tpA, xh, j2rp2, j4rp4)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(:, :), INTENT(IN) :: xh
               TYPE(helio_pl), INTENT(INOUT)              :: helio_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb(npl, helio_plA, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: npl
               REAL(DP), INTENT(IN)           :: dt
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_kickvb
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb_tp(ntp, helio_tpA, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: ntp
               REAL(DP), INTENT(IN)           :: dt
               TYPE(helio_tp), INTENT(INOUT)  :: helio_tpA
          END SUBROUTINE helio_kickvb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift(npl, swiftest_plA, dt, pt)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl
               REAL(DP), INTENT(IN)                   :: dt
               REAL(DP), DIMENSION(:), INTENT(OUT) :: pt
               TYPE(swiftest_pl), INTENT(INOUT)       :: swiftest_plA
          END SUBROUTINE helio_lindrift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift_tp(ntp, swiftest_tpA, dt, pt)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), INTENT(IN)                  :: dt
               REAL(DP), DIMENSION(:), INTENT(IN) :: pt
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE helio_lindrift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step(lfirst, lextra_force, t, npl, ntp, helio_plA, helio_tpA, j2rp2, j4rp4, dt)
               USE swiftest_globals
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)      :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)   :: lfirst
               INTEGER(I4B), INTENT(IN)      :: npl, ntp
               REAL(DP), INTENT(IN)          :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE helio_step
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_pl(lfirst, lextra_force, t, npl, helio_plA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                 :: lfirst
               INTEGER(I4B), INTENT(IN)                    :: npl
               REAL(DP), INTENT(IN)                        :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(:), INTENT(OUT)         :: ptb, pte
               REAL(DP), DIMENSION(:,:), INTENT(OUT)       :: xbeg,xend
               TYPE(helio_pl), INTENT(INOUT)               :: helio_plA
          END SUBROUTINE helio_step_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_tp(lfirsttp, lextra_force, t, npl, ntp, helio_plA, helio_tpA, j2rp2, j4rp4, dt, &
               xbeg, xend, ptb, pte)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                :: lfirsttp
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(:), INTENT(IN)         :: ptb, pte
               REAL(DP), DIMENSION(:,:), INTENT(IN)       :: xbeg, xend
               TYPE(helio_pl), INTENT(INOUT)              :: helio_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_step_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch(t, npl, helio_plA)
               USE swiftest_globals
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
          END SUBROUTINE helio_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch_tp(t, ntp, helio_tpA)
               USE swiftest_globals
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: ntp
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE helio_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergesub, symba_plA, & 
               discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                       :: lbig_discard
               INTEGER(I4B), INTENT(IN)                       :: npl, nsppl, nsptp, nmergesub
               REAL(DP), INTENT(IN)                           :: t, mtiny
               CHARACTER(*), INTENT(IN)                       :: fname
               TYPE(symba_pl), INTENT(INOUT)                  :: symba_plA
               TYPE(swiftest_tp), INTENT(INOUT)               :: discard_tpA
               TYPE(swiftest_pl), INTENT(INOUT)               :: discard_plA
               TYPE(symba_merger), INTENT(INOUT)              :: mergeadd_list, mergesub_list
          END SUBROUTINE io_discard_write_symba
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_pl(npl, swiftest_plA, param)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: npl
               TYPE(swiftest_pl), INTENT(INOUT):: swiftest_plA
               type(user_input_parameters),intent(inout) :: param

          END SUBROUTINE io_dump_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_tp(ntp, swiftest_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: ntp
               TYPE(swiftest_tp), INTENT(INOUT):: swiftest_tpA
          END SUBROUTINE io_dump_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_tp(intpfile, in_type, ntp, symba_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               CHARACTER(*), INTENT(IN)         :: intpfile, in_type
               TYPE(symba_tp), INTENT(INOUT)    :: symba_tpA
          END SUBROUTINE io_init_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_open(iu, fname, fopenstat, fmt, ierr)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: iu
               INTEGER(I4B), INTENT(OUT) :: ierr
               CHARACTER(*), INTENT(IN)  :: fname, fopenstat, fmt
          END SUBROUTINE io_open
     END INTERFACE

     INTERFACE
          FUNCTION io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B)                           :: io_read_encounter
               INTEGER(I4B), INTENT(OUT)              :: name1, name2
               REAL(DP), INTENT(OUT)                  :: t, mass1, mass2
               REAL(DP), DIMENSION(:), INTENT(OUT)    :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)               :: encounter_file
          END FUNCTION io_read_encounter
     END INTERFACE

     INTERFACE
          FUNCTION io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B)               :: io_read_hdr
               INTEGER(I4B), INTENT(IN)   :: iu
               INTEGER(I4B), INTENT(OUT)  :: npl, ntp, iout_form
               REAL(DP), INTENT(OUT)      :: t
               CHARACTER(*), INTENT(IN)   :: out_type
          END FUNCTION io_read_hdr
     END INTERFACE

     INTERFACE
          FUNCTION io_read_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B)                    :: io_read_line
               INTEGER(I4B), INTENT(IN)        :: iu
               INTEGER(I4B), INTENT(OUT)       :: name
               REAL(DP), INTENT(OUT)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(OUT) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)        :: out_type
          END FUNCTION io_read_line
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
               xh1, xh2, vh1, vh2, encounter_file)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: name1, name2
               REAL(DP), INTENT(IN)                  :: t, mass1, mass2, radius1, radius2
               REAL(DP), DIMENSION(:), INTENT(IN)    :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)              :: encounter_file
          END SUBROUTINE io_write_encounter
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: iu, name
               REAL(DP), INTENT(IN)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(IN) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)       :: out_type
          END SUBROUTINE io_write_line
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_acc(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                    :: npl
               REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4
               REAL(DP), DIMENSION(:), INTENT(IN)        :: irh
               REAL(DP), DIMENSION(:,:), INTENT(IN)  :: xh
               REAL(DP), DIMENSION(:,:), INTENT(OUT) :: aobl
               TYPE(swiftest_pl), INTENT(INOUT)            :: swiftest_plA
          END SUBROUTINE obl_acc
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                    :: ntp
               REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4, msun
               REAL(DP), DIMENSION(:), INTENT(IN)        :: irht
               REAL(DP), DIMENSION(:, :), INTENT(IN)  :: xht
               REAL(DP), DIMENSION(:, :), INTENT(OUT) :: aoblt
          END SUBROUTINE obl_acc_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_pot(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, oblpot)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl
               REAL(DP), INTENT(IN)                       :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                      :: oblpot
               REAL(DP), DIMENSION(:), INTENT(IN)       :: irh
               REAL(DP), DIMENSION(:, :), INTENT(IN) :: xh
               TYPE(swiftest_pl), INTENT(INOUT)           :: swiftest_plA
          END SUBROUTINE obl_pot
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_scget(angle, sx, cx)
            !!$omp declare simd(orbel_scget)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: angle
               REAL(DP), INTENT(OUT) :: sx, cx
          END SUBROUTINE orbel_scget
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aeq(x, v, mu, a, e, q)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(:), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, q
          END SUBROUTINE orbel_xv2aeq
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(:), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, q, capm, tperi
          END SUBROUTINE orbel_xv2aqt
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(:), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, inc, capom, omega, capm
          END SUBROUTINE orbel_xv2el
     END INTERFACE

     INTERFACE
          SUBROUTINE python_io_write_frame_pl(t, symba_plA, npl, out_stat)
               use swiftest_globals
               use swiftest_data_structures
               use module_helio
               use module_symba
               IMPLICIT NONE
               real(DP), intent(in)      :: t
               type(symba_pl),intent(in) :: symba_plA
               integer, intent(in)       :: npl
               character(*), intent(in)  :: out_stat
          END SUBROUTINE python_io_write_frame_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)
               use swiftest_globals
               use swiftest_data_structures
               use module_helio
               use module_symba
               IMPLICIT NONE
               real(DP), intent(in)      :: t
               type(symba_tp),intent(in) :: symba_tpA
               integer, intent(in)       :: ntp
               character(*), intent(in)  :: out_stat
          END SUBROUTINE python_io_write_frame_tp
     END INTERFACE



     INTERFACE
          SUBROUTINE rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(:), INTENT(IN) :: xr, vr
               INTEGER(I4B), INTENT(OUT)             :: iflag
          END SUBROUTINE rmvs_chk_ind
     END INTERFACE

     INTERFACE
      function symba_casedisruption (symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, &
                                     mass_res, param, Qloss) result(status)
         use swiftest_globals
         use swiftest_data_structures
         use module_symba
         implicit none
         type(symba_pl), intent(inout)             :: symba_plA
         integer(I4B), dimension(:), intent(in)    :: idx_parent
         integer(I4B), intent(inout)               :: nmergeadd
         type(symba_merger), intent(inout)         :: mergeadd_list
         real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
         real(DP), dimension(:),   intent(in)      :: mass, radius, mass_res
         type(user_input_parameters),intent(inout) :: param
         real(DP), intent(in)                      :: Qloss
         integer(I4B)                              :: status
      end function symba_casedisruption

      function symba_casehitandrun (symba_plA, idx_parent, nmergeadd, mergeadd_list, name, x, v, mass, radius, Lspin, Ip, &
                                    mass_res, param, Qloss) result(status)
         use swiftest_globals
         use swiftest_data_structures
         use module_symba
         implicit none
         type(symba_pl), intent(inout)             :: symba_plA
         integer(I4B), dimension(:), intent(in)    :: idx_parent
         integer(I4B), intent(inout)               :: nmergeadd
         type(symba_merger), intent(inout)         :: mergeadd_list
         integer(I4B), dimension(:), intent(in)    :: name
         real(DP), dimension(:,:), intent(in)      :: x, v, Lspin, Ip
         real(DP), dimension(:), intent(in)        :: mass, radius, mass_res
         type(user_input_parameters),intent(inout) :: param
         real(DP), intent(in)                      :: Qloss
         integer(I4B)                              :: status
      end function symba_casehitandrun

      function symba_casemerge (symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, lspin, Ip, param) result(status)
         use swiftest_globals
         use swiftest_data_structures
         use module_symba
         implicit none
         type(symba_pl), intent(inout)             :: symba_plA
         integer(I4B), dimension(:), intent(in)    :: idx_parent
         integer(I4B), intent(inout)               :: nmergeadd
         type(symba_merger), intent(inout)         :: mergeadd_list
         real(DP), dimension(:,:), intent(in)      :: x, v, lspin, Ip
         real(DP), dimension(:), intent(in)        :: mass, radius
         type(user_input_parameters),intent(inout) :: param
         integer(I4B)                              :: status
      end function symba_casemerge

      function symba_casesupercatastrophic (symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, lspin, Ip, &
                                                mass_res, param, Qloss) result(status)
         use swiftest_globals
         use swiftest_data_structures
         use module_symba
         implicit none
         type(symba_pl), intent(inout)             :: symba_plA
         integer(I4B), dimension(:), intent(in)    :: idx_parent
         integer(I4B), intent(inout)               :: nmergeadd
         type(symba_merger), intent(inout)         :: mergeadd_list
         real(DP), dimension(:,:), intent(in)      :: x, v, lspin, Ip
         real(DP), dimension(:), intent(in)        :: mass, radius, mass_res
         type(user_input_parameters),intent(inout) :: param
         real(DP), intent(in)                      :: Qloss
         integer(I4B)                              :: status
      end function symba_casesupercatastrophic
     end interface

     INTERFACE
          SUBROUTINE symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(OUT)          :: lencounter, lvdotr
               INTEGER(I4B), INTENT(IN)           :: irec
               REAL(DP), INTENT(IN)               :: rhill1, rhill2, dt
               REAL(DP), DIMENSION(:), INTENT(IN) :: xr, vr
          END SUBROUTINE symba_chk
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_collision(t, symba_plA, nplplenc, plplenc_list, ldiscard, mergeadd_list, nmergeadd, param)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               real(DP), intent(in)                      :: t
               integer(I4B), intent(inout)               :: nplplenc, nmergeadd
               type(symba_pl)                            :: symba_plA
               type(symba_plplenc), intent(inout)        :: plplenc_list
               type(symba_merger), intent(inout)         :: mergeadd_list
               logical, intent(inout)                    :: ldiscard
               type(user_input_parameters),intent(inout) :: param
          END SUBROUTINE symba_collision
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_pl(t, dt, index_enc, nmergesub, mergesub_list, npl, symba_plA, plplenc_list, param)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                  :: index_enc
               INTEGER(I4B), INTENT(IN)                  :: npl
               INTEGER(I4B), INTENT(INOUT)               :: nmergesub
               REAL(DP), INTENT(IN)                      :: t, dt
               TYPE(symba_plplenc), INTENT(INOUT)        :: plplenc_list
               TYPE(symba_merger), INTENT(INOUT)         :: mergesub_list
               TYPE(symba_pl), INTENT(INOUT)             :: symba_plA
               TYPE(user_input_parameters),intent(inout) :: param
          END SUBROUTINE symba_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscard)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl
               REAL(DP), INTENT(IN)          :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)      :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               LOGICAL(LGT), INTENT(INOUT)   :: ldiscard
          END SUBROUTINE symba_discard_peri_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_pl(t, npl, ntp, symba_plA, symba_tpA, rmin, rmax, rmaxu, qmin, qmin_coord,          &
               qmin_alo, qmin_ahi, ldiscard)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT)    :: npl, ntp
               REAL(DP), INTENT(IN)           :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)       :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)  :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)  :: symba_tpA
               LOGICAL(LGT), INTENT(INOUT)    :: ldiscard
          END SUBROUTINE symba_discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_sun_pl(t, npl, ntp, msys, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, ldiscard)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl, ntp
               REAL(DP), INTENT(IN)             :: t, msys, rmin, rmax, rmaxu
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               LOGICAL(LGT), INTENT(INOUT)      :: ldiscard
          END SUBROUTINE symba_discard_sun_pl
     END INTERFACE

     INTERFACE
            subroutine symba_discard_conserve_mtm(swiftest_plA, ipl, lescape)
            use swiftest_globals
            use swiftest_data_structures
            implicit none
            integer(I4B), intent(in)    :: ipl
            type(swiftest_pl), intent(inout) :: swiftest_plA
            logical, intent(in)        :: lescape
            end subroutine
      END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_tp(t, npl, ntp, symba_plA, symba_tpA, dt, &
               rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, lrhill_present, ldiscard_tp)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)      :: lrhill_present
               INTEGER(I4B), INTENT(IN)      :: npl
               INTEGER(I4B), INTENT(INOUT)   :: ntp
               REAL(DP), INTENT(IN)          :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)      :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               TYPE(symba_tp), INTENT(INOUT) :: symba_tpA
               LOGICAL(LGT), INTENT(INOUT)   :: ldiscard_tp
          END SUBROUTINE symba_discard_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_energy(npl, symba_plA, j2rp2, j4rp4, ke_orbit, ke_spin, pe, te, Ltot)
               USE swiftest_globals
               use module_symba
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)            :: npl
               REAL(DP), INTENT(IN)                :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)               :: ke_orbit, ke_spin, pe, te
               REAL(DP), DIMENSION(:), INTENT(OUT) :: Ltot
               TYPE(symba_pl), INTENT(INOUT)       :: symba_plA
          END SUBROUTINE symba_energy
     END INTERFACE

     INTERFACE
         subroutine symba_frag_pos (param, symba_plA, idx_parent, x, v, L_spin, Ip, mass, radius, &
                                    Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, lmerge, Qloss)
            use swiftest_globals
            USE swiftest_data_structures
            USE module_symba
            implicit none
            type(user_input_parameters), intent(in)   :: param 
            type(symba_pl), intent(inout)             :: symba_plA
            integer(I4B), dimension(:), intent(in)    :: idx_parent
            real(DP), intent(in)                      :: Qloss
            real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
            real(DP), dimension(:), intent(in)        :: mass, radius, m_frag, rad_frag
            real(DP), dimension(:,:), intent(in)      :: Ip_frag
            real(DP), dimension(:,:), intent(out)     :: xb_frag, vb_frag, rot_frag
            logical, intent(out)                      :: lmerge
         end subroutine symba_frag_pos
      END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, symba_plA, j2rp2, j4rp4, nplplenc, &
               plplenc_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplplenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch_tp(lextra_force, t, npl, nplm, ntp, symba_plA, symba_tpA, &
               xh, j2rp2, j4rp4,  &
               npltpenc, pltpenc_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, ntp, npltpenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(:, :), INTENT(IN)    :: xh
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
               TYPE(symba_pltpenc), INTENT(IN)               :: pltpenc_list
          END SUBROUTINE symba_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift(irec, npl, symba_plA, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: irec, npl
               REAL(DP), INTENT(IN)          :: dt
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
          END SUBROUTINE symba_helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift_tp(irec, ntp, symba_tpA, mu, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: irec, ntp
               REAL(DP), INTENT(IN)          :: mu, dt
               TYPE(symba_tp), INTENT(INOUT) :: symba_tpA
          END SUBROUTINE symba_helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch(lflag, lextra_force, t, npl, nplm, helio_plA, j2rp2, j4rp4)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)       :: npl, nplm
               REAL(DP), INTENT(IN)           :: t, j2rp2, j4rp4
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE symba_helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch_int(npl, nplm, helio_plA)
               USE swiftest_globals
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl, nplm
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
          END SUBROUTINE symba_helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_plA, &
               symba_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: irec, nplplenc, npltpenc
               REAL(DP), INTENT(IN)            :: dt, sgn
               TYPE(symba_plplenc), INTENT(IN) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(IN) :: pltpenc_list
               TYPE(symba_pl), INTENT(INOUT)   :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)   :: symba_tpA
          END SUBROUTINE symba_kick
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_tp(t, dt, index_enc, pltpenc_list, vbs, encounter_file, symba_plA, symba_tpA)
               USE swiftest_globals
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index_enc
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file
               TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          END SUBROUTINE symba_merge_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lfirst
               INTEGER(I4B), INTENT(IN)       :: npl
               REAL(DP), INTENT(IN)           :: msys
               CHARACTER(*), INTENT(IN)       :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)  :: symba_plA
          END SUBROUTINE symba_peri
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_rearray(npl, nplm, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
                                   discard_tpA, ldiscard, ldiscard_tp, mtiny)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT)                   :: npl, nplm, ntp, nsppl, nsptp, nmergeadd 
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
               TYPE(symba_tp), INTENT(INOUT)              :: discard_tpA
               TYPE(symba_pl), INTENT(INOUT)              :: discard_plA
               TYPE(symba_merger), INTENT(INOUT)             :: mergeadd_list 
               LOGICAL(LGT), INTENT(IN)                      :: ldiscard, ldiscard_tp 
               real(DP), intent(in)                          :: mtiny
          END SUBROUTINE symba_rearray

     END INTERFACE  


     INTERFACE
          SUBROUTINE symba_reorder_pl(npl, symba_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(symba_pl), INTENT(INOUT)  :: symba_plA
               INTEGER(I4B)                              :: i
               INTEGER(I4B), DIMENSION(:), ALLOCATABLE   :: index
               REAL(DP), DIMENSION(:), ALLOCATABLE       :: mass
               REAL(DP), DIMENSION(:,:), allocatable     :: symba_plwkspA
               INTEGER(I4B), DIMENSION(:,:), allocatable :: symba_plwkspA_id_status
          END SUBROUTINE symba_reorder_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swiftest_pl1P, &
               swiftest_tp1P)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                            :: npl, ntp
               TYPE(swiftest_pl), POINTER                           :: swiftest_pl1P
               TYPE(swiftest_tp), POINTER                           :: swiftest_tp1P
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               TYPE(symba_tp), INTENT(INOUT) :: symba_tpA
               TYPE(symba_pl), POINTER                             :: symba_pl1P
               TYPE(symba_tp), POINTER                             :: symba_tp1P
          END SUBROUTINE symba_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step(t, dt, param,npl, ntp,symba_plA, symba_tpA,  &
                         nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, &
                         mergeadd_list, mergesub_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               use user
               IMPLICIT NONE
               TYPE(user_input_parameters), INTENT(INOUT)       :: param        ! Derived type containing user defined parameters 
               INTEGER(I4B), INTENT(IN)                         :: npl, ntp
               INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, dt
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio(lfirst, lextra_force, t, npl, nplm, ntp, helio_plA, helio_tpA, j2rp2, &
               j4rp4, dt)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)      :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)   :: lfirst
               INTEGER(I4B), INTENT(IN)      :: npl, nplm, ntp
               REAL(DP), INTENT(IN)          :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE symba_step_helio
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, helio_plA, j2rp2, j4rp4, dt, xbeg, xend,    &
               ptb, pte)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                     :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                  :: lfirst
               INTEGER(I4B), INTENT(IN)                     :: npl, nplm
               REAL(DP), INTENT(IN)                         :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(:,:), INTENT(OUT)        :: xbeg, xend
               REAL(DP), DIMENSION(:), INTENT(OUT)          :: ptb, pte
               TYPE(helio_pl), INTENT(INOUT)                :: helio_plA
          END SUBROUTINE symba_step_helio_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_interp(t, npl, nplm, ntp, symba_plA, symba_tpA,  &
               dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,    &
               mergesub_list, param)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               use user
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)           :: npl, nplm, ntp, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, dt
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
               type(user_input_parameters), intent(inout) :: param
          END SUBROUTINE symba_step_interp
     END INTERFACE

     INTERFACE
          RECURSIVE SUBROUTINE symba_step_recur(t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, &
            nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, param)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               use user
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)           :: ireci, npl, nplm, ntp, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, dt0
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
               type(user_input_parameters), intent(inout)  :: param
          END SUBROUTINE symba_step_recur
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch(t, npl, symba_plA)
               USE swiftest_globals
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)     :: npl
               REAL(DP), INTENT(IN)         :: t
               TYPE(symba_pl), INTENT(INOUT):: symba_plA
          END SUBROUTINE symba_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch_tp(t, ntp, symba_tpA)
               USE swiftest_globals
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: ntp
               REAL(DP), INTENT(IN)           :: t
               TYPE(symba_tp), INTENT(INOUT)  :: symba_tpA
          END SUBROUTINE symba_user_getacch_tp
     END INTERFACE

     INTERFACE
         SUBROUTINE symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke_orbit, ke_spin, pe, te, Ltot)
               USE swiftest_globals
               USE swiftest_data_structures
               use module_symba
               IMPLICIT NONE
               integer(I4B), intent(in)              :: npl
               type(symba_pl), intent(inout)         :: symba_plA
               real(DP), intent(in)                  :: j2rp2, j4rp4
               real(DP), intent(out)                 :: ke_orbit, ke_spin, pe, te
               real(DP), dimension(:), intent(out)   :: Ltot
         END SUBROUTINE symba_energy_eucl
     END INTERFACE

     INTERFACE 
          subroutine symba_chk_eucl(npl, irec, symba_plA, dt, plplenc_list, nplplenc)
               use swiftest_globals
               use swiftest_data_structures
               use module_symba
               IMPLICIT NONE
               integer(I4B), intent(in)                         :: npl, irec
               type(symba_pl), intent(inout)                    :: symba_plA
               real(DP), intent(in)                             :: dt
               type(symba_plplenc), intent(inout)  :: plplenc_list
               integer(I4B), intent(inout)                      :: nplplenc
          end subroutine symba_chk_eucl
     END INTERFACE

     INTERFACE 
          SUBROUTINE symba_chk_eucl_pltp(symba_plA, symba_tpA, dt, lencounter, lvdotr, npltpenc)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               IMPLICIT NONE
               TYPE(symba_pl), INTENT(IN)                    :: symba_plA
               TYPE(symba_tp), INTENT(IN)                    :: symba_tpA
               REAL(DP), INTENT(IN)               :: dt
               LOGICAL(LGT), DIMENSION(:), INTENT(OUT) :: lencounter, lvdotr
               INTEGER(I4B), INTENT(INOUT)        :: npltpenc
          END SUBROUTINE symba_chk_eucl_pltp
     END INTERFACE


     INTERFACE
          SUBROUTINE symba_getacch_eucl(lextra_force, t, npl, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplplenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_getacch_eucl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch_tp_eucl(lextra_force, t, npl, ntp, symba_plA, symba_tpA, &
               xh, j2rp2, j4rp4, npltpenc, pltpenc_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, ntp, npltpenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(:, :), INTENT(IN)         :: xh
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
               TYPE(symba_pltpenc), INTENT(IN)               :: pltpenc_list
          END SUBROUTINE symba_getacch_tp_eucl
     END INTERFACE


     INTERFACE
          SUBROUTINE symba_step_eucl(t,dt,param,npl, ntp,symba_plA, symba_tpA,  &
            nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, &
            mergeadd_list, mergesub_list)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_helio
               USE module_symba
               use user
               IMPLICIT NONE
               TYPE(user_input_parameters), INTENT(INOUT)       :: param        ! Derived type containing user defined parameters 
               INTEGER(I4B), INTENT(IN)                         :: npl, ntp
               INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, dt
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_eucl
     END INTERFACE


     INTERFACE
         SUBROUTINE symba_step_interp_eucl(t, npl, nplm, ntp, symba_plA, symba_tpA,&
            dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, &
            mergeadd_list, mergesub_list, param)
               USE swiftest_globals
               USE swiftest_data_structures
               USE module_symba
               USE user
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, nplm, ntp, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)                :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                       :: t, dt
               TYPE(symba_pl), INTENT(INOUT)              :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)              :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT)         :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT)         :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)          :: mergeadd_list, mergesub_list
               type(user_input_parameters), intent(inout) :: param
          END SUBROUTINE symba_step_interp_eucl
     END INTERFACE

     interface
         subroutine symba_step_reset(npl, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list)
         use swiftest_globals
         use swiftest_data_structures
         USE module_symba
         implicit none
         integer(I4B), intent(in) :: npl
         type(symba_pl), intent(inout)  :: symba_plA
         type(symba_tp), intent(inout)  :: symba_tpA
         type(symba_plplenc), intent(inout)   :: plplenc_list
         type(symba_pltpenc), intent(inout)   :: pltpenc_list
         type(symba_merger), intent(inout)    :: mergeadd_list, mergesub_list
         end subroutine symba_step_reset
      end interface 

     interface
         function symba_mergeadd_eoffset(npl, symba_plA, mergeadd_list, mergesub_list, addi, &
            addf, subi, subf, param) result(eoffset)
            USE swiftest_globals
            USE swiftest_data_structures
            USE module_symba
            USE user
            implicit none
            integer(I4B), intent(in)                :: npl
            type(symba_pl), intent(in)           :: symba_plA
            type(symba_merger), intent(in)          :: mergeadd_list, mergesub_list
            integer(I4B), intent(in)                :: addi, addf, subi, subf
            type(user_input_parameters), intent(in) :: param
            real(DP)                                :: eoffset
         end function symba_mergeadd_eoffset
      end interface


      INTERFACE
         SUBROUTINE util_dist_index_plpl(npl, nplm, symba_plA)
            USE swiftest_globals
            USE swiftest_data_structures
            use module_symba
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)  :: npl, nplm
            type(symba_pl), intent(inout) :: symba_plA
         END SUBROUTINE
      END INTERFACE

   INTERFACE
     SUBROUTINE util_dist_index_pltp(nplm, ntp, symba_tpA)
          USE swiftest_globals
          USE swiftest_data_structures
          use module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)  :: nplm, ntp
          type(symba_tp), intent(inout) :: symba_tpA
     END SUBROUTINE util_dist_index_pltp
   END INTERFACE

INTERFACE
     SUBROUTINE util_dist_eucl_plpl(invar, outvar, symba_plA)
          USE swiftest_globals
          USE swiftest_data_structures
          use module_symba
          IMPLICIT NONE
          REAL(DP),DIMENSION(:,:),INTENT(IN) :: invar
          REAL(DP), DIMENSION(:,:),INTENT(INOUT) :: outvar
          type(symba_pl), intent(inout) :: symba_plA
     END SUBROUTINE util_dist_eucl_plpl
END INTERFACE

INTERFACE
     SUBROUTINE util_dist_eucl_pltp(planets, test_particles, outvar, symba_tpA)
          USE swiftest_globals
          USE swiftest_data_structures
          use module_symba
          IMPLICIT NONE
          REAL(DP),DIMENSION(:,:),INTENT(IN) :: planets
          REAL(DP),DIMENSION(:,:),INTENT(IN) :: test_particles
          REAL(DP), DIMENSION(:,:),INTENT(INOUT) :: outvar
          type(symba_tp), intent(inout) :: symba_tpA
     END SUBROUTINE
END INTERFACE


     INTERFACE 
          pure SUBROUTINE util_crossproduct(ar1, ar2, ans)
               USE swiftest_globals
               IMPLICIT NONE
               real(DP),dimension(:),intent(in)  :: ar1,ar2
               real(DP),dimension(:),intent(out) :: ans
          END SUBROUTINE
     END INTERFACE

     INTERFACE
          SUBROUTINE util_exit(code)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: code
          END SUBROUTINE util_exit
     END INTERFACE

     INTERFACE
          SUBROUTINE util_hills(npl, swiftest_plA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE util_hills
     END INTERFACE

     INTERFACE
          SUBROUTINE util_index(arr, index)
               USE swiftest_globals
               USE module_nrutil
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
               REAL(DP), DIMENSION(:), INTENT(IN)      :: arr
          END SUBROUTINE util_index
     END INTERFACE

     INTERFACE
          SUBROUTINE util_peri(lfirst, ntp, swiftest_tpA, mu, msys, qmin_coord)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)         :: lfirst
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, msys
               CHARACTER(*), INTENT(IN)         :: qmin_coord
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE util_peri
     END INTERFACE

     INTERFACE 
          SUBROUTINE util_resize_pl(symba_plA, npl_new, npl_old)
               USE swiftest_globals
               USE module_symba
               USE swiftest_data_structures
               USE module_helio
               USE module_nrutil
               USE module_swiftestalloc
               IMPLICIT NONE
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               INTEGER(I4B), INTENT(IN)      :: npl_old, npl_new
          END SUBROUTINE util_resize_pl
     END INTERFACE

     INTERFACE util_sort
          SUBROUTINE util_sort_i4b(arr)
               USE swiftest_globals
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_i4b
          SUBROUTINE util_sort_sp(arr)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_sp
          SUBROUTINE util_sort_dp(arr)
               USE swiftest_globals
               IMPLICIT NONE
               REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_dp
     END INTERFACE

     INTERFACE
          SUBROUTINE util_toupper(string)
               USE swiftest_globals
               IMPLICIT NONE
               CHARACTER(*), INTENT(INOUT) :: string
          END SUBROUTINE util_toupper
     END INTERFACE

     INTERFACE
          SUBROUTINE util_valid(npl, ntp, swiftest_plA, swiftest_tpA)
               USE swiftest_globals
               USE swiftest_data_structures
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl, ntp
               TYPE(swiftest_pl), INTENT(IN) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(IN) :: swiftest_tpA
          END SUBROUTINE util_valid
     END INTERFACE

     INTERFACE
          SUBROUTINE util_version
               USE swiftest_globals
               IMPLICIT NONE
          END SUBROUTINE util_version
     END INTERFACE

     ! Added by D. Minton
     INTERFACE
         FUNCTION util_kahan_sum(xsum_current, xi, xerror) 
            USE swiftest_globals
            IMPLICIT NONE
            REAL(DP)                :: util_kahan_sum
            REAL(DP), INTENT(IN)    :: xsum_current, xi
            REAL(DP), INTENT(INOUT) :: xerror
         END FUNCTION
     END INTERFACE

     INTERFACE
         FUNCTION collresolve_resolve(model,m1,m2,r1,r2,p1,p2,v1,v2,n,mres,rres,pres,vres)
         USE swiftest_globals
         IMPLICIT NONE
         INTEGER(I4B) :: collresolve_resolve
         INTEGER(I4B), INTENT(IN) :: model               ! collision model to apply
         REAL(DP), INTENT(IN) :: m1                      ! mass of the target
         REAL(DP), INTENT(IN) :: m2                      ! mass of the impactor
         REAL(DP), INTENT(IN) :: r1                      ! radius of the target
         REAL(DP), INTENT(IN) :: r2                      ! radius of the impactor
         REAL(DP), DIMENSION(:), INTENT(IN) :: p1        ! position of the target
         REAL(DP), DIMENSION(:), INTENT(IN) :: p2        ! position of the impactor
         REAL(DP), DIMENSION(:), INTENT(IN) :: v1        ! velocity of the target
         REAL(DP), DIMENSION(:), INTENT(IN) :: v2        ! velocity of the impactor
         INTEGER, INTENT(IN) :: n                                ! number of bodies to return
         REAL(DP), DIMENSION(:), INTENT(OUT) :: mres   ! mass of the resulting bodies
         REAL(DP), DIMENSION(:), INTENT(OUT) :: rres   ! radius of the resulting bodies
         REAL(DP), DIMENSION(:,:), INTENT(OUT) :: pres ! position of the resulting bodies
         REAL(DP), DIMENSION(:,:), INTENT(OUT) :: vres ! velocity of the resulting bodies
         END FUNCTION
      END INTERFACE

     INTERFACE
         SUBROUTINE symba_regime(Mcb, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, mtiny, Qloss)
          USE swiftest_globals
          USE module_symba
          USE swiftest_data_structures
          USE module_helio
          USE module_nrutil
          USE module_swiftestalloc
         IMPLICIT NONE
          INTEGER(I4B), INTENT(OUT)              :: regime
          REAL(DP), INTENT(INOUT)                :: Mcb, Mlr, Mslr, m1, m2, rad1, rad2, den1, den2, mtiny
          REAL(DP), DIMENSION(:), INTENT(IN)     :: xh1, xh2, vb1, vb2
          real(DP), intent(out)                  :: Qloss
         END SUBROUTINE symba_regime
     END INTERFACE

END MODULE module_interfaces
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


    
