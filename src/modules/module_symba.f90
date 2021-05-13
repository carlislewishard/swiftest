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

   type, public, extends(lambda_obj) :: symba_vel_lambda_obj
      procedure(abstract_objective_func), pointer, nopass :: ke_objective_func_ptr => null()
      real(DP), dimension(:),   allocatable :: m_frag, v_r_mag
      real(DP), dimension(:,:), allocatable :: v_r_unit
      real(DP), dimension(NDIM)             :: L_lin_tan
      real(DP)                              :: T_rad
   contains
      generic   :: init => ke_objective_func_init
      procedure :: eval => ke_objective_func_eval
      procedure :: ke_objective_func_init
      final     :: ke_objective_func_destroy
   end type symba_vel_lambda_obj

   abstract interface
      function abstract_objective_func(v_r_mag_unknowns, v_r_mag_knowns, m_frag, v_r_unit, L_lin_tan, T_rad) result(fnorm)
         ! Template for the kinetic energy constraint function used for minimizing
         import DP
         real(DP), dimension(:),   intent(in) :: v_r_mag_unknowns   !! Unknown radial velocity magnitudes
         real(DP), dimension(:),   intent(in) :: v_r_mag_knowns  !! Known Radial velocity magnitude
         real(DP), dimension(:),   intent(in) :: m_frag    !! Fragment masses
         real(DP), dimension(:,:), intent(in) :: v_r_unit  !! Radial unit vectors
         real(DP), dimension(:),   intent(in) :: L_lin_tan !! Tangential component of linear momentum
         real(DP),                 intent(in) :: T_rad     !! Target radial kinetic energ
         real(DP)                             :: fnorm     !! The objective function result: norm of the vector composed of the tangential momentum and energy
      end function
   end interface

   contains
      subroutine ke_objective_func_init(self, lambda, v_r_mag, m_frag, v_r_unit, L_lin_tan, T_rad)
         implicit none
         ! Arguments
         class(symba_vel_lambda_obj), intent(out) :: self
         procedure(abstract_objective_func)       :: lambda  !! The lambda function
         real(DP), dimension(:),      intent(in)  :: v_r_mag   !! Radial velocity magnitude
         real(DP), dimension(:),      intent(in)  :: m_frag    !! Fragment masses
         real(DP), dimension(:,:),    intent(in)  :: v_r_unit  !! Radial unit vectors
         real(DP), dimension(:),      intent(in)  :: L_lin_tan !! Tangential component of linear momentum
         real(DP),                    intent(in)  :: T_rad     !! Target radial kinetic ener
   
         self%ke_objective_func_ptr  => lambda
         allocate(self%m_frag, source=m_frag)
         allocate(self%v_r_mag, source=v_r_mag)
         allocate(self%v_r_unit, source=v_r_unit)
         self%L_lin_tan(:) = L_lin_tan(:)
         self%T_rad = T_rad
      end subroutine ke_objective_func_init

      subroutine ke_objective_func_destroy(self)
         implicit none
         type(symba_vel_lambda_obj) :: self
         if (allocated(self%m_frag)) deallocate(self%m_frag)
         if (allocated(self%v_r_unit)) deallocate(self%v_r_unit)
         if (allocated(self%v_r_mag)) deallocate(self%v_r_mag)
         if (associated(self%ke_objective_func_ptr)) nullify(self%ke_objective_func_ptr)
      end subroutine ke_objective_func_destroy 

   function ke_objective_func_eval(self, x) result(fnorm)
      implicit none
      ! Arguments
      class(symba_vel_lambda_obj), intent(in) :: self
      real(DP), dimension(:),      intent(in) :: x
      ! Result
      real(DP)                      :: fnorm
      ! Internals
      integer(I4B) :: nfrag, nunknown

      if (associated(self%ke_objective_func_ptr)) then
         nfrag = size(self%v_r_mag)
         nunknown = size(x)
         fnorm = self%ke_objective_func_ptr(x, self%v_r_mag(nunknown + 1:nfrag), self%m_frag, self%v_r_unit, self%L_lin_tan, self%T_rad)
      else
         error stop "KE Objective function was not initialized."
      end if
   end function ke_objective_func_eval
END MODULE module_symba