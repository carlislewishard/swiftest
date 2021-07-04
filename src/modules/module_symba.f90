!**********************************************************************************************************************************
!
!  Unit Name   : module_symba
!  Unit Type   : module
!  Project   : SWIFTEST
!  Package   : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Symplectic Massive Body Algorithm
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
MODULE module_symba

   USE swiftest_globals
   USE module_helio
   use lambda_function
   IMPLICIT NONE

   INTEGER(I4B), PARAMETER :: NENMAX = 32767
   INTEGER(I4B), PARAMETER :: NTENC = 3
   REAL(DP), PARAMETER   :: RHSCALE = 6.5_DP
   REAL(DP), PARAMETER   :: RSHELL = 0.48075_DP

   type symba_kinship
      integer(I4B) :: parent ! Index of parent particle
      integer(I4B) :: nchild ! number of children in merger list
      integer(I4B), dimension(:), allocatable :: child ! Index of children particles
   end type symba_kinship

   type symba_pl
      logical(LGT), dimension(:),   allocatable :: lcollision ! flag indicating whether body has merged with another this time step
      logical(LGT), dimension(:),   allocatable :: lencounter! flag indicating whether body is part of an encounter this time step
      integer(I4B), dimension(:),   allocatable :: nplenc  ! number of encounters with other planets this time step
      integer(I4B), dimension(:),   allocatable :: ntpenc  ! number of encounters with test particles this time step
      integer(I4B), dimension(:),   allocatable :: levelg  ! level at which this body should be moved
      integer(I4B), dimension(:),   allocatable :: levelm  ! deepest encounter level achieved this time step
      integer(I4B), dimension(:),   allocatable :: isperi  ! perihelion passage flag
      real(DP),   dimension(:),     allocatable :: peri    ! perihelion distance
      real(DP),   dimension(:),     allocatable :: atp   ! semimajor axis following perihelion passage
      type(helio_pl)                            :: helio   ! HELIO planet structure
      type(symba_kinship), dimension(:), allocatable :: kin  ! Array of merger relationship structures that can account for multiple pairwise 
                                                             ! mergers in a single step
   end type symba_pl

   type symba_tp
      integer(I4B), dimension(:),   allocatable :: nplenc  ! number of encounters with planets this time step
      integer(I4B), dimension(:),   allocatable :: levelg  ! level at which this particle should be moved
      integer(I4B), dimension(:),   allocatable :: levelm  ! deepest encounter level achieved this time step
      type(helio_tp)                  :: helio   ! HELIO test particle structure
   end type symba_tp

   type symba_plplenc
      integer(I4B)                              :: nplplenc ! Total number of pl-pl encounters
      logical(LGT), dimension(:),   allocatable :: lvdotr ! relative vdotr flag
      integer(I4B), dimension(:),   allocatable :: status ! status of the interaction
      integer(I4B), dimension(:),   allocatable :: level  ! encounter recursion level
      integer(I4B), dimension(:),   allocatable :: index1   ! position of the first planet in encounter
      integer(I4B), dimension(:),   allocatable :: index2   ! position of the second planet in encounter
      integer(I4B), dimension(:),   allocatable :: enc_child   ! the child of the encounter
      integer(I4B), dimension(:),   allocatable :: enc_parent   ! the child of the encounter
      real(DP),     dimension(:,:), allocatable :: xh1          ! the heliocentric position of parent 1 in encounter
      real(DP),     dimension(:,:), allocatable :: xh2          ! the heliocentric position of parent 2 in encounter
      real(DP),     dimension(:,:), allocatable :: vb1          ! the barycentric velocity of parent 1 in encounter
      real(DP),     dimension(:,:), allocatable :: vb2          ! the barycentric velocity of parent 2 in encounter
   end type symba_plplenc

   type symba_pltpenc
      integer(I4B)                              :: npltpenc ! Total number of pl-tp encounters
      logical(LGT), dimension(:),   allocatable :: lvdotr ! relative vdotr flag
      integer(I4B), dimension(:),   allocatable :: status ! status of the interaction
      integer(I4B), dimension(:),   allocatable :: level  ! encounter recursion level
      integer(I4B), dimension(:),   allocatable :: indexpl    ! position of the planet in encounter
      integer(I4B), dimension(:),   allocatable :: indextp    ! position of the test particle in encounter
   end type symba_pltpenc

   type symba_merger
      integer(I4B), dimension(:),   allocatable :: id   ! external identifier
      integer(I4B), dimension(:),   allocatable :: index_ps ! position of the particle
      integer(I4B), dimension(:),   allocatable :: status   ! status
      integer(I4B), dimension(:),   allocatable :: nadded   ! number of resultant bodies from this collisional event aka 
                                        !      number of fragments
      real(DP),   dimension(:,:),   allocatable :: xb     ! barycentric position
      real(DP),   dimension(:,:),   allocatable :: vb     ! barycentric velocity
      real(DP),   dimension(:),     allocatable :: mass   ! mass
      real(DP),   dimension(:),     allocatable :: radius   ! radius
      real(DP),   dimension(:,:),   allocatable :: IP     ! moment of intertia
      real(DP),   dimension(:,:),   allocatable :: rot    ! rotation
      type(swiftest_particle_info), dimension(:), allocatable :: info

   end type symba_merger 

   interface
      module subroutine symba_read_pl_in(symba_plA, param) 
         type(symba_pl),              intent(inout) :: symba_plA  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_parameters), intent(inout) :: param    !! Input collection of user-defined parameters
      end subroutine symba_read_pl_in
   end interface

END MODULE module_symba