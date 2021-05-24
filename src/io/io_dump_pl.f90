!**********************************************************************************************************************************
!
!  Unit Name   : io_dump_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Dump planet data to file
!
!  Input
!    Arguments : npl            : number of planets
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                lclose         : logical flag indicating whether to check for planet-test particle encounters
!                lrhill_present : logical flag indicating whether Hill's sphere radii are present in planet data
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : to dump file
!                npl            : number of planets
!                id             : planet identifier     (from planet structure, for each planet)
!                mass           : mass                  (from planet structure, for each planet)
!                rhill          : Hill's sphere radius  (from planet structure, for each planet except the Sun, if lrhill_present)
!                radius         : planet radius         (from planet structure, for each planet except the Sun, if lclose)
!                xh             : heliocentric position (from planet structure, for each planet)
!                vh             : heliocentric velocity (from planet structure, for each planet)
!
!  Invocation  : CALL io_dump_pl(npl, swifter_pl1P, lclose, lrhill_present)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_dump_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE io_dump_pl(npl, swiftest_plA, param)

! Modules
     USE swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => io_dump_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)         :: npl
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
      type(user_input_parameters),intent(inout) :: param

! Internals
   INTEGER(I4B)                     :: ierr
   INTEGER(I4B), SAVE               :: idx = 1
   integer(I4B),parameter             :: LUN = 7

   open(unit = LUN, file = DUMP_PL_FILE(idx), form = "UNFORMATTED", status = 'REPLACE', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Unable to open binary dump file ", trim(dump_pl_file(idx))
      call util_exit(failure)
   end if
   write(LUN) npl
   write(LUN) swiftest_plA%id(1:npl)
   write(LUN) swiftest_plA%mass(1:npl)
   if (param%lrhill_present) write(LUN) swiftest_plA%rhill(1:npl) 
   if (param%lclose) write(LUN) swiftest_plA%radius(1:npl) 
   write(LUN) swiftest_plA%xh(:,1:npl)
   write(LUN) swiftest_plA%vh(:,1:npl)
   if (param%lrotation) THEN
      write(LUN) swiftest_plA%Ip(:,1:npl)
      write(LUN) swiftest_plA%rot(:,1:npl)
   end if
   close(LUN)
   idx = idx + 1
   if (idx > 2) idx = 1

   return

end subroutine io_dump_pl
!**********************************************************************************************************************************
!
!  author(s)   : david e. kaufmann
!
!  revision control system (rcs) information
!
!  source file : $rcsfile$
!  full path   : $source$
!  revision    : $revision$
!  date        : $date$
!  programmer  : $author$
!  locked by   : $locker$
!  state       : $state$
!
!  modification history:
!
!  $log$
!**********************************************************************************************************************************
