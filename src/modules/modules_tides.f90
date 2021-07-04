!**********************************************************************************************************************************
!
!  Unit Name   : module_tides
!  Unit Type   : module
!  Project   : SWIFTEST
!  Package   : module
!  Language    : Fortran 90/95
!  author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
!  Description : Definition of data and structures specific to the calculation of tidal torques
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File    : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File    : N/A
!
!  Invocation  : N/A
!
!  Notes     : 
!
!**********************************************************************************************************************************
MODULE module_tides

   USE swiftest_globals
   USE module_helio
   USE module_symba
   use lambda_function
   IMPLICIT NONE

   INTEGER(I4B), PARAMETER :: NENMAX = 32767
   INTEGER(I4B), PARAMETER :: NTENC = 3
   REAL(DP), PARAMETER     :: RHSCALE = 6.5_DP
   REAL(DP), PARAMETER     :: RSHELL = 0.48075_DP

   type tides
      !integer(I4B), dimension(:),   allocatable :: id   ! external identifier

      !type(swiftest_particle_info), dimension(:), allocatable :: info

   end type tides 

   interface
      module subroutine symba_read_pl_in(symba_plA, param) 
         type(symba_pl),              intent(inout) :: symba_plA  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_parameters), intent(inout) :: param    !! Input collection of user-defined parameters
      end subroutine symba_read_pl_in
   end interface

END MODULE module_tides