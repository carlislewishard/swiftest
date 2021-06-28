module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   implicit none

   type swiftest_particle_info
      sequence
      character(len=32)         :: origin_type 
      real(DP)                  :: origin_time 
      real(DP), dimension(NDIM) :: origin_xh
      real(DP), dimension(NDIM) :: origin_vh
   end type swiftest_particle_info

   type,public :: swiftest_tp
      integer(I4B)                                :: nbody = 0  !! Number of bodies
      integer(I4B)                                :: maxid = 0  !! Maximum id
      character(len=STRMAX), dimension(:),  allocatable :: name   !! Non-unique name
      integer(I4B), dimension(:),     allocatable :: id     !! External identifier (unique)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      integer(I4B), dimension(:),     allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),     allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),     allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
      integer(I4B), dimension(:,:),   allocatable :: k_pltp
      integer(I8B)                                :: num_pltp_comparisons
      type(swiftest_particle_info), dimension(:), allocatable :: info
   contains
      procedure :: alloc => swiftest_tp_allocate
      procedure :: dealloc => swiftest_tp_deallocate
      procedure :: read_from_file => swiftest_read_tp_in 
   end type swiftest_tp

   type,public,extends(swiftest_tp) :: swiftest_pl
      real(DP),     dimension(:),     allocatable :: mass   !! Mass
      real(DP),     dimension(:),     allocatable :: radius !! Radius
      real(DP),     dimension(:),     allocatable :: rhill  !! Hill's sphere radius
      real(DP),     dimension(:,:),   allocatable :: Ip     !! Unitless principal moments of inertia (I1, I2, I3) / (MR**2). Principal axis rotation assumed. 
      real(DP),     dimension(:,:),   allocatable :: rot    !! Body rotation vector in inertial coordinate frame 
      real(DP),     dimension(NDIM)               :: Lcb_initial  !! Initial angular momentum of the central body
      real(DP),     dimension(NDIM)               :: dLcb = [0.0_DP, 0.0_DP, 0.0_DP] !! Change in angular momentum of the central body
      real(DP),     dimension(NDIM)               :: Lescape = [0.0_DP, 0.0_DP, 0.0_DP] !! Angular momentum of bodies that escaped the system (used for bookeeping)
      real(DP)                                    :: Mcb_initial !! Initial mass of the central body
      real(DP)                                    :: dMcb = 0.0_DP !! Change in mass of the central body
      real(DP)                                    :: Mescape = 0.0_DP !! Mass of bodies that escaped the system (used for bookeeping)
      real(DP)                                    :: Rcb_initial !! Initial radius of the central body
      real(DP)                                    :: dRcb = 0.0_DP!! Change in the radius of the central body
      real(DP)                                    :: Ecollisions = 0.0_DP !! Energy lost from system due to collisions
      real(DP)                                    :: Euntracked = 0.0_DP !! Energy gained from system due to escaped bodies
      integer(I4B), dimension(:,:), allocatable   :: k_plpl
      integer(I8B)                                :: num_plpl_comparisons
   contains
      procedure :: alloc => swiftest_pl_allocate
      procedure :: dealloc => swiftest_pl_deallocate
      procedure :: read_from_file => swiftest_read_pl_in 
   end type swiftest_pl

   type, public :: swiftest_parameters
      real(DP)             :: t0 = 0.0_DP          !! Integration start time
      real(DP)             :: tstop = 0.0_DP       !! Integration stop time
      real(DP)             :: dt = 0.0_DP          !! Time step
      character(STRMAX)    :: inplfile = ''        !! Name of input file for planets
      character(STRMAX)    :: intpfile = ''        !! Name of input file for test particles
      character(STRMAX)    :: in_type = 'ASCII'    !! Format of input data files
      integer(I4B)         :: istep_out = -1       !! Number of time steps between binary outputs
      character(STRMAX)    :: outfile = ''         !! Name of output binary file
      character(STRMAX)    :: particle_file = ''        !! Name of output particle information file
      character(STRMAX)    :: out_type = REAL4_TYPE!! Binary format of output file
      character(STRMAX)    :: out_form = 'XV'      !! Data to write to output file
      character(STRMAX)    :: out_stat = 'NEW'     !! Open status for output binary file
      integer(I4B)         :: istep_dump = -1      !! Number of time steps between dumps
      real(DP)             :: j2rp2 = 0.0_DP       !! J2 * R**2 for the Sun
      real(DP)             :: j4rp4 = 0.0_DP       !! J4 * R**4 for the Sun
      real(DP)             :: rmin = -1.0_DP       !! Minimum heliocentric radius for test particle
      real(DP)             :: rmax = -1.0_DP       !! Maximum heliocentric radius for test particle
      real(DP)             :: rmaxu = -1.0_DP      !! Maximum unbound heliocentric radius for test particle
      real(DP)             :: qmin = -1.0_DP       !! Minimum pericenter distance for test particle
      character(STRMAX)    :: qmin_coord = 'HELIO' !! Coordinate frame to use for qmin
      real(DP)             :: qmin_alo = -1.0_DP   !! Minimum semimajor axis for qmin
      real(DP)             :: qmin_ahi = -1.0_DP   !! Maximum semimajor axis for qmin
      character(STRMAX)    :: encounter_file = ''  !! Name of output file for encounters
      real(DP)             :: mtiny = 0.0_DP       !! Smallest mass that is fully gravitating
      character(STRMAX)    :: ring_outfile = ''    !! Name of output file in ring moons
      real(DP)             :: MU2KG = -1.0_DP      !! Converts mass units to grams
      real(DP)             :: TU2S  = -1.0_DP      !! Converts time units to seconds
      real(DP)             :: DU2M = -1.0_DP       !! Converts distance unit to centimeters
      integer(I4B), dimension(:), allocatable :: seed  !! Random seeds


      !Logical flags to turn on or off various features of the code
      logical :: lextra_force = .false.            !! User defined force function turned on
      logical :: lbig_discard = .false.            !! Save big bodies on every discard
      logical :: lrhill_present = .false.          !! Hill's radius is in input file
      logical :: lclose = .false.                  !! Turn on close encounters
      logical :: lfragmentation = .false.          !! Do fragmentation modeling instead of simple merger.
      logical :: lmtiny     = .false.              !! Use the MTINY variable (Automatically set if running SyMBA)
      logical :: lrotation  = .false.              !! Include rotation states of big bodies
      logical :: ltides     = .false.              !! Include tidal dissipation 
      logical :: lringmoons = .false.              !! Turn on the ringmoons code 
      logical :: lenergy = .false.                 !! Track the total energy of the system
      logical :: lfirstenergy = .true.
      real(DP)                 :: Eorbit_orig = 0.0_DP
      real(DP)                 :: Mtot_orig = 0.0_DP
      real(DP)                 :: Lmag_orig = 0.0_DP
      real(DP), dimension(NDIM) :: Ltot_orig = 0.0_DP
      real(DP), dimension(NDIM) :: Lorbit_orig = 0.0_DP
      real(DP), dimension(NDIM) :: Lspin_orig = 0.0_DP
      logical :: lfirstkick = .true.               !! Initiate the first kick in a symplectic step

      ! Future features not implemented or in development
      logical :: lgr = .false.               !! Turn on GR
      logical :: lyarkovsky = .false.        !! Turn on Yarkovsky effect
      logical :: lyorp = .false.             !! Turn on YORP effect
   contains
      procedure :: read_from_file => io_read_param_in
      procedure :: dump_to_file => io_dump_param
      procedure :: udio_reader => io_udio_reader
      procedure :: udio_writer => io_udio_writer
      !TODO: Figure out if user-defined derived-type io can be made to work properly
      !generic   :: read(formatted) => udio_reader
      !generic   :: write(formatted) => udio_writer
   end type swiftest_parameters

   interface

      module subroutine swiftest_read_pl_in(self, param) 
         class(swiftest_pl),          intent(inout) :: self  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_parameters), intent(inout) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_pl_in

      module subroutine swiftest_read_tp_in(self, param) 
         class(swiftest_tp),          intent(inout) :: self  !! Swiftest data structure to store massive body initial conditions
         type(swiftest_parameters), intent(inout) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_tp_in

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function io_get_token

      !> Interface for type-bound procedure to read in the input parameters from a file
      module subroutine io_read_param_in(param, inparfile, swiftest_plA) 
         class(swiftest_parameters),intent(out) :: param         !! Input collection of user-defined parameters
         character(*), intent(in)                 :: inparfile     !! Parameter input file name (i.e. param.in)
         type(swiftest_pl), intent(inout)         :: swiftest_plA
      end subroutine io_read_param_in

      !> Interface for type-bound procedure to write out the user parameters into a dump file in case the run needs to be restarted
      module subroutine io_dump_param(param, t, swiftest_plA)
         class(swiftest_parameters),intent(in)  :: param    !! Output collection of user-defined parameters
         real(DP),intent(in)                      :: t        !! Current simulation time
         type(swiftest_pl), intent(inout)         :: swiftest_plA
      end subroutine io_dump_param

      !> Interface for type-bound procedure for user-defined derived-type IO for reading
      module subroutine io_udio_reader(param, unit, iotype, v_list, iostat, iomsg, swiftest_plA) 
         class(swiftest_parameters),intent(inout)  :: param         !! Input collection of user-defined parameters
         integer, intent(in)                    :: unit        
         character(len=*), intent(in)           :: iotype
         integer, intent(in)                    :: v_list(:)
         integer, intent(out)                   :: iostat
         character(len=*), intent(inout)        :: iomsg
         type(swiftest_pl), intent(inout)       :: swiftest_plA
      end subroutine io_udio_reader

      !> Interface for type-bound procedure for user-defined derived-type IO for writing
      module subroutine io_udio_writer(param, unit, iotype, v_list, iostat, iomsg, swiftest_plA) 
         class(swiftest_parameters),intent(in)  :: param         !! Output collection of user-defined parameters
         integer, intent(in)                 :: unit        
         character(len=*), intent(in)        :: iotype
         integer, intent(in)                 :: v_list(:)
         integer, intent(out)                :: iostat
         character(len=*), intent(inout)     :: iomsg
         type(swiftest_pl), intent(inout)    :: swiftest_plA
      end subroutine io_udio_writer

   end interface

   contains

      subroutine swiftest_tp_allocate(self,n)
         implicit none

         class(swiftest_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate

         self%nbody = n
         if (n <= 0) return
         allocate(self%id(n))
         allocate(self%status(n))
         allocate(self%peri(n))
         allocate(self%atp(n))
         allocate(self%isperi(n))
         allocate(self%xh(NDIM,n))
         allocate(self%vh(NDIM,n))
         allocate(self%xb(NDIM,n))
         allocate(self%vb(NDIM,n))
         allocate(self%info(n))

         self%id = 0
         self%status = 0
         self%peri = 0.0_DP
         self%atp = 0.0_DP
         self%isperi = 0.0_DP
         self%xh = 0.0_DP
         self%vh = 0.0_DP
         self%xb = 0.0_DP
         self%vb = 0.0_DP
         return
      end subroutine swiftest_tp_allocate

      subroutine swiftest_pl_allocate(self,n)
         implicit none

         class(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
         integer, intent(in)                  :: n    !! Number of massive bodies to allocate

         self%nbody = n
         if (n <= 0) return
         call self%swiftest_tp%alloc(n)

         allocate(self%mass(n))
         allocate(self%radius(n))
         allocate(self%rhill(n))
         allocate(self%Ip(NDIM,n))
         allocate(self%rot(NDIM,n))

         self%mass = 0.0_DP
         self%radius = 0.0_DP
         self%rhill = 0.0_DP
         self%Ip = 0.0_DP
         self%rot = 0.0_DP
         return
      end subroutine swiftest_pl_allocate

      subroutine swiftest_tp_deallocate(self)
         implicit none

         class(swiftest_tp), intent(inout)    :: self

         self%nbody = 0
         if (allocated(self%id)) deallocate(self%id)
         if (allocated(self%status)) deallocate(self%status)
         if (allocated(self%isperi)) deallocate(self%isperi)
         if (allocated(self%peri)) deallocate(self%peri)
         if (allocated(self%atp)) deallocate(self%atp)
         if (allocated(self%xh)) deallocate(self%xh)
         if (allocated(self%vh)) deallocate(self%vh)
         if (allocated(self%xb)) deallocate(self%xb)
         if (allocated(self%vb)) deallocate(self%vb)
         if (allocated(self%k_pltp)) deallocate(self%k_pltp)
         if (allocated(self%info)) deallocate(self%info)
         return
      end subroutine swiftest_tp_deallocate

      subroutine swiftest_pl_deallocate(self)
         implicit none

         class(swiftest_pl), intent(inout)    :: self

         call self%swiftest_tp%dealloc()
         if (allocated(self%mass)) deallocate(self%mass)
         if (allocated(self%radius)) deallocate(self%radius)
         if (allocated(self%rhill)) deallocate(self%rhill)
         if (allocated(self%Ip)) deallocate(self%Ip)
         if (allocated(self%rot)) deallocate(self%rot)
         if (allocated(self%k_plpl)) deallocate(self%k_plpl)
         return
      end subroutine swiftest_pl_deallocate


end module swiftest_data_structures
