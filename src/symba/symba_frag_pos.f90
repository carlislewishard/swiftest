module symba_frag_pos_lambda_implementation
   use swiftest_globals
   use lambda_function
   implicit none

   ! Define the class and interface used to implement the lambda function
   type, public, extends(lambda_obj) :: ke_constraint
      procedure(abstract_objective_func), pointer, nopass :: ke_objective_func_ptr => null()
      real(DP), dimension(:),   allocatable :: v_t_mag
      real(DP), dimension(:,:), allocatable :: v_r_unit, v_t_unit
      real(DP), dimension(:),   allocatable :: m_frag
      real(DP), dimension(NDIM)             :: vcom
      real(DP)                              :: ke_target
   contains
      generic   :: init => ke_objective_func_init
      procedure :: eval => ke_objective_func_eval
      procedure, nopass :: ke_objective_func_init
      final     :: ke_objective_func_destroy
   end type ke_constraint
   interface ke_constraint
      module procedure ke_objective_func_init
   end interface

   abstract interface
      function abstract_objective_func(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom, ke_target) result(fnorm)
         ! Template for the kinetic energy constraint function used for minimizing
         import DP
         real(DP), dimension(:),   intent(in) :: v_r_mag   !! Radial velocity mangitude
         real(DP), dimension(:,:), intent(in) :: v_r_unit  !! Radial velocity unit vector
         real(DP), dimension(:),   intent(in) :: v_t_mag   !! Tangential velocity magnitude
         real(DP), dimension(:,:), intent(in) :: v_t_unit  !! Tangential velocity unit vector
         real(DP), dimension(:),   intent(in) :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in) :: vcom      !! Barycentric velocity of collisional system center of mass
         real(DP),                 intent(in) :: ke_target !! Target kinetic energ
         real(DP)                             :: fnorm     !! The objective function result: norm of the vector composed of the tangential momentum and energy
      end function
   end interface

   contains
      type(ke_constraint) function ke_objective_func_init(lambda, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom, ke_target)
         implicit none
         ! Arguments
         procedure(abstract_objective_func)   :: lambda    !! The lambda function
         real(DP), dimension(:,:), intent(in) :: v_r_unit  !! Radial velocity unit vector
         real(DP), dimension(:),   intent(in) :: v_t_mag   !! Tangential velocity magnitude
         real(DP), dimension(:,:), intent(in) :: v_t_unit  !! Tangential velocity unit vector
         real(DP), dimension(:),   intent(in) :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in) :: vcom      !! Barycentric velocity of collisional system center of mass
         real(DP),                 intent(in) :: ke_target !! Target kinetic energ
         ! Internals
         associate(self => ke_objective_func_init)
            self%ke_objective_func_ptr  => lambda
            allocate(self%v_r_unit, source=v_r_unit)
            allocate(self%v_t_mag, source=v_t_mag)
            allocate(self%v_t_unit, source=v_t_unit)
            allocate(self%m_frag, source=m_frag)
            self%vcom(:) = vcom(:)
            self%ke_target = ke_target
         end associate
         return
      end function ke_objective_func_init

      subroutine ke_objective_func_destroy(self)
         implicit none
         type(ke_constraint) :: self
         if (allocated(self%v_r_unit)) deallocate(self%v_r_unit)
         if (allocated(self%v_t_mag)) deallocate(self%v_t_mag)
         if (allocated(self%v_t_unit)) deallocate(self%v_t_unit)
         if (allocated(self%m_frag)) deallocate(self%m_frag)
         if (associated(self%ke_objective_func_ptr)) nullify(self%ke_objective_func_ptr)
      end subroutine ke_objective_func_destroy 

   function ke_objective_func_eval(self, x) result(fval) 
      implicit none
      ! Arguments
      class(ke_constraint),   intent(in) :: self
      real(DP), dimension(:), intent(in) :: x
      ! Result
      real(DP)                           :: fval

      if (associated(self%ke_objective_func_ptr)) then
         fval = self%ke_objective_func_ptr(x, self%v_r_unit, self%v_t_mag, self%v_t_unit, self%m_frag, self%vcom, self%ke_target)
      else
         error stop "KE Objective function was not initialized."
      end if
   end function ke_objective_func_eval

end module symba_frag_pos_lambda_implementation

subroutine symba_frag_pos(param, symba_plA, family, x, v, L_spin, Ip, mass, radius, &
                          nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lmerge)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Places the collision fragments on a circle oriented with a plane defined
   !! by the position and velocity vectors of the collision
   !! 
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, except_this_one => symba_frag_pos
   implicit none
   ! Arguments
   type(user_input_parameters), intent(in) :: param 
   type(symba_pl), intent(inout)           :: symba_plA
   integer(I4B), dimension(:), intent(in)  :: family
   real(DP), dimension(:,:), intent(inout) :: x, v, L_spin, Ip
   real(DP), dimension(:),   intent(inout) :: mass, radius
   integer(I4B), intent(inout)             :: nfrag
   real(DP), dimension(:),   intent(inout) :: m_frag, rad_frag
   real(DP), dimension(:,:), intent(inout) :: Ip_frag
   real(DP), dimension(:,:), intent(inout) :: xb_frag, vb_frag, rot_frag
   logical, intent(out)                    :: lmerge ! Answers the question: Should this have been a merger instead?
   real(DP), intent(inout)                 :: Qloss
   ! Internals
   real(DP), parameter                     :: f_spin = 0.20_DP !! Fraction of pre-impact orbital angular momentum that is converted to fragment spin
   real(DP)                                :: mscale = 1.0_DP, rscale = 1.0_DP, vscale = 1.0_DP, tscale = 1.0_DP, Lscale = 1.0_DP, Escale = 1.0_DP! Scale factors that reduce quantities to O(~1) in the collisional system
   real(DP)                                :: mtot 
   real(DP), dimension(NDIM)               :: xcom, vcom
   integer(I4B)                            :: ii, fam_size
   logical, dimension(:), allocatable      :: lexclude
   real(DP), dimension(NDIM, 2)            :: rot, L_orb 
   real(DP), dimension(:,:), allocatable   :: x_frag, v_frag, v_r_unit, v_t_unit
   real(DP), dimension(:), allocatable     :: rmag, v_r_mag, v_t_mag
   real(DP), dimension(NDIM)               :: Ltot_before
   real(DP), dimension(NDIM)               :: Ltot_after
   real(DP)                                :: Etot_before, ke_orbit_before, ke_spin_before, pe_before, Lmag_before
   real(DP)                                :: Etot_after,  ke_orbit_after,  ke_spin_after,  pe_after,  Lmag_after, dEtot, dLmag
   real(DP), dimension(NDIM)               :: L_frag_spin, L_frag_tot, L_frag_orb, L_offset
   real(DP)                                :: ke_family, ke_target, ke_frag, ke_offset
   real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
   character(len=*), parameter             :: fmtlabel = "(A14,10(ES11.4,1X,:))"
   integer(I4B)                            :: try, ntry
   integer(I4B), parameter                 :: NFRAG_MIN = 6 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
   real(DP)                                :: r_max_start = 1.0_DP
   real(DP), parameter                     :: Ltol = 2 * epsilon(1.0_DP)
   real(DP), parameter                     :: Etol = 1e-8_DP
   integer(I4B), parameter                 :: MAXTRY = 1000
   integer(I4B) :: iflip = 1

   if (nfrag < NFRAG_MIN) then
      write(*,*) "symba_frag_pos needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
      lmerge = .true.
      return
   end if

   ntry = nfrag - NFRAG_MIN

   allocate(x_frag, source=xb_frag)
   allocate(v_frag, source=vb_frag)

   associate(npl => symba_plA%helio%swiftest%nbody, status => symba_plA%helio%swiftest%status) 
      allocate(lexclude(npl))
      where (status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
         lexclude(1:npl) = .true.  
      elsewhere
         lexclude(1:npl) = .false. 
      end where
   end associate
   call set_scale_factors()

   call calculate_system_energy(linclude_fragments=.false.)


   !write(*,        "(' -------------------------------------------------------------------------------------')")
   !write(*,        "('              Energy normalized by |Etot_before|')")
   
   call define_coordinate_system()
   call define_pre_collisional_family()
  
   try = 1
   lmerge = .false.
   do while (nfrag >= NFRAG_MIN .and. try <= MAXTRY)
      lmerge = .false.
      !write(*,        "(' -------------------------------------------------------------------------------------')")
      !write(*,*) "Try number: ",try, ' of ',ntry
      !write(*,        "(' -------------------------------------------------------------------------------------')")
      !write(*,        "('  First pass to get angular momentum ')")
   
      call set_fragment_position_vectors()
      call set_fragment_tangential_velocities()
      call calculate_system_energy(linclude_fragments=.true.)
   
      !write(*,        "(' -------------------------------------------------------------------------------------')")
      !write(*,        "('             |      T_orb      T_spin           T          pe        Etot        Ltot')")
      !write(*,        "(' -------------------------------------------------------------------------------------')")
      !write(*,fmtlabel) ' change      |',(ke_orbit_after - ke_orbit_before) / abs(Etot_before), &
      !                                    (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
      !                                    (ke_orbit_after + ke_spin_after - ke_orbit_before - ke_spin_before)/ abs(Etot_before), &
      !                                    (pe_after - pe_before) / abs(Etot_before), &
      !                                    dEtot, dLmag
  ! 
  !    write(*,        "(' -------------------------------------------------------------------------------------')")
  !    write(*,        "(' If T_offset > 0, this will be a merge')") 
  !    write(*,fmtlabel) ' T_offset  |',ke_offset / abs(Etot_before)
   
      call set_fragment_radial_velocities(lmerge)
   
      call calculate_system_energy(linclude_fragments=.true.)
   
      if (.not.lmerge) then
         if (dEtot > 0.0_DP)  then
            write(*,*) 'Failed try ',try,': Positive energy'
            lmerge = .true.
         else if (abs((Etot_after - Etot_before + Qloss) / Etot_before) > Etol) then
            write(*,*) 'Failed try ',try,': Energy error too big: ',dEtot
            lmerge = .true.
         else if (abs(dLmag) > Ltol) then
            write(*,*) 'Failed try ',try,': Angular momentum error too big: ',dLmag
            lmerge = .true.
         end if

         !write(*,        "(' -------------------------------------------------------------------------------------')")
         !write(*,        "('  Second pass to get energy ')")
         !write(*,        "(' -------------------------------------------------------------------------------------')")
         !write(*,fmtlabel) ' T_frag targ |',ke_target / abs(Etot_before)
         !write(*,fmtlabel) ' T_frag calc |',ke_frag / abs(Etot_before)
         !write(*,        "(' T_offset should now be small')") 
         !write(*,fmtlabel) ' T_offset    |',ke_offset / abs(Etot_before)
         !write(*,        "(' -------------------------------------------------------------------------------------')")
         !write(*,        "('             |      T_orb      T_spin           T          pe        Etot        Ltot')")
         !write(*,        "(' -------------------------------------------------------------------------------------')")
         !write(*,fmtlabel) ' change      |',(ke_orbit_after - ke_orbit_before) / abs(Etot_before), &
         !                                    (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
         !                                    (ke_orbit_after + ke_spin_after - ke_orbit_before - ke_spin_before)/ abs(Etot_before), &
         !                                    (pe_after - pe_before) / abs(Etot_before), &
         !                                    dEtot, dLmag
      !else
         !write(*,*) 'Failed try ',try,': Could not find radial velocities.'
      end if
 
      if (.not.lmerge) exit
      call restructure_failed_fragments()
      try = try + 1
   end do
   write(*,        "(' -------------------------------------------------------------------------------------')")
   write(*,        "('  Final diagnostic')")
   write(*,        "(' -------------------------------------------------------------------------------------')")
   if (lmerge) then
      write(*,*) "symba_frag_pos failed after: ",try," tries"
      stop
   else
      write(*,*) "symba_frag_pos succeeded after: ",try," tries"
      write(*,        "(' dL_tot should be very small' )")
      write(*,fmtlabel) ' dL_tot      |', dLmag
      write(*,        "(' dE_tot should be negative and equal to Qloss' )")
      write(*,fmtlabel) ' dE_tot      |', dEtot
      write(*,fmtlabel) ' Qloss       |', -Qloss / abs(Etot_before)
      write(*,fmtlabel) ' dE - Qloss  |', (Etot_after - Etot_before + Qloss) / abs(Etot_before)
   end if
   write(*,        "(' -------------------------------------------------------------------------------------')")

   call restore_scale_factors()

   return 

   contains

   ! Because of the complexity of this procedure, we have chosen to break it up into a series of nested subroutines.

   subroutine set_scale_factors()
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. This scaling makes it easier for the non-linear minimization 
      !! to converge on a solution
      implicit none

      ! Find the center of mass of the collisional system	
      mtot = sum(mass(:)) 
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

      ! Set scale factors
      mscale = mtot !! Because of the implied G, mass is actually G*mass with units of distance**3 / time**2
      rscale = sum(radius(:))
      vscale = sqrt(mscale / rscale)
      tscale = rscale / vscale
      Lscale = mscale * rscale * vscale
      Escale = mscale * vscale**2

      xcom(:) = xcom(:) / rscale
      vcom(:) = vcom(:) / vscale

      mtot = mtot / mscale
      mass = mass / mscale
      radius = radius / rscale
      x = x / rscale
      v = v / vscale
      L_spin = L_spin / Lscale

      m_frag = m_frag / mscale
      rad_frag = rad_frag / rscale
      Qloss = Qloss / Escale

      return
   end subroutine set_scale_factors

   subroutine restore_scale_factors()
      !! author: David A. Minton
      !!
      !! Restores dimenional quantities back to the system units
      implicit none
      integer(I4B) :: i

      mtot = mscale
      mass = mass * mscale
      radius = radius * rscale
      x = x * rscale
      v = v * vscale
      L_spin = L_spin * Lscale

      x_frag = x_frag * rscale
      v_frag = v_frag * vscale
      m_frag = m_frag * mscale
      rot_frag = rot_frag / tscale
      rad_frag = rad_frag * rscale
      ! Convert the final result to the system units
      xcom(:) = xcom(:) * rscale
      vcom(:) = vcom(:) * vscale
      do i = 1, nfrag
         xb_frag(:, i) = x_frag(:, i) + xcom(:)
         vb_frag(:, i) = v_frag(:, i) + vcom(:)
      end do

      Qloss = Qloss * Escale
      mscale = 1.0_DP
      rscale = 1.0_DP
      vscale = 1.0_DP
      tscale = 1.0_DP
      Lscale = 1.0_DP
      Escale = 1.0_DP

      return
   end subroutine restore_scale_factors

   subroutine define_coordinate_system()
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate system, including the unit vectors of both the system and individual fragments.
      implicit none
      integer(I4B) :: i
      real(DP), dimension(NDIM) ::  x_cross_v, xc, vc, delta_r, delta_v
      real(DP)   :: r_col_norm, v_col_norm

      allocate(rmag(nfrag))
      allocate(v_r_mag(nfrag))
      allocate(v_t_mag(nfrag))
      allocate(v_r_unit(NDIM,nfrag))
      allocate(v_t_unit(NDIM,nfrag))

      L_orb(:, :) = 0.0_DP
      ! Compute orbital angular momentum of pre-impact system
      do i = 1, 2
         xc(:) = x(:, i) - xcom(:) 
         vc(:) = v(:, i) - vcom(:)
         call util_crossproduct(xc(:), vc(:), x_cross_v(:))
         L_orb(:, i) = mass(i) * x_cross_v(:)
      end do

      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      L_frag_tot(:) = L_spin(:,1) + L_spin(:,2) + L_orb(:, 1) + L_orb(:, 2)
      L_frag_spin(:) = L_spin(:,1) + L_spin(:, 2) + f_spin * (L_orb(:, 1) + L_orb(:, 2))
      L_frag_orb(:) =  L_frag_tot - L_frag_spin

      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))     
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))

      ! We will initialize fragments on a plane defined by the pre-impact system, with the z-axis aligned with the angular momentum vector
      ! and the y-axis aligned with the pre-impact distance vector.
      y_col_unit(:) = delta_r(:) / r_col_norm  
      z_col_unit(:) = L_frag_tot(:) / norm2(L_frag_tot)
      ! The cross product of the z- by x-axis will give us the y-axis
      call util_crossproduct(y_col_unit, z_col_unit, x_col_unit)

      return
   end subroutine define_coordinate_system

   subroutine define_pre_collisional_family()
      !! author: David A. Minton
      !!
      !! Defines the pre-collisional "family" consisting of all of the bodies involved in the collision. These bodies need to be identified so that they can be excluded from energy calculations once
      !! the collisional products (the fragments) are created.
      integer(I4B) :: i
      associate(vbpl => symba_plA%helio%swiftest%vb, Mpl => symba_plA%helio%swiftest%mass)
         fam_size = size(family)

         ! We need the original kinetic energy of just the pre-impact family members in order to balance the energy later
         ke_family = 0.0_DP
         do i = 1, fam_size
            ke_family = ke_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i)))
            lexclude(family(i)) = .true. ! For all subsequent energy calculations the pre-impact family members will be replaced by the fragments
         end do
         ke_family = 0.5_DP * ke_family / Escale
      end associate
      return
   end subroutine define_pre_collisional_family

   subroutine calculate_system_energy(linclude_fragments)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy, including all bodies in the symba_plA list that do not have a corresponding value of the lexclude array that is true
      !! and optionally including fragments.
      use module_swiftestalloc
      implicit none
      ! Arguments
      logical,                intent(in) :: linclude_fragments
      ! Internals
      real(DP) :: ke_orbit, ke_spin, pe, te
      real(DP), dimension(NDIM)  :: Lorbit, Lspin
      integer(I4B) :: i, npl_new, nplm
      logical, dimension(:), allocatable :: ltmp
      logical :: lk_plpl
      type(symba_pl) :: symba_plwksp

      ! Because we're making a copy of symba_pl with the excludes/fragments appended, we need to deallocate the
      ! big k_plpl array and recreate it when we're done, otherwise we run the risk of blowing up the memory by
      ! allocating two of these ginormous arrays simulteouously. This is not particularly efficient, but as this
      ! subroutine should be called relatively infrequently, it shouldn't matter too much.
      !if (allocated(symba_plA%helio%swiftest%k_plpl)) deallocate(symba_plA%helio%swiftest%k_plpl) 

      ! Build the internal planet list out of the non-excluded bodies and optionally with fragments appended. This
      ! will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it
      ! is in the main program.
      associate(npl => symba_plA%helio%swiftest%nbody)
         lk_plpl = allocated(symba_plA%helio%swiftest%k_plpl)
         if (lk_plpl) deallocate(symba_plA%helio%swiftest%k_plpl) 
         if (linclude_fragments) then ! Temporarily expand the planet list to feed it into symba_energy
            npl_new = npl + nfrag
         else
            npl_new  = npl
         end if
         call symba_pl_allocate(symba_plwksp, npl_new)
   
         ! Copy over old data
         symba_plwksp%helio%swiftest%id(1:npl) = symba_plA%helio%swiftest%id(1:npl)
         symba_plwksp%helio%swiftest%status(1:npl) = symba_plA%helio%swiftest%status(1:npl)
         symba_plwksp%helio%swiftest%mass(1:npl) = symba_plA%helio%swiftest%mass(1:npl) / mscale
         symba_plwksp%helio%swiftest%radius(1:npl) = symba_plA%helio%swiftest%radius(1:npl) / rscale
         symba_plwksp%helio%swiftest%xh(:,1:npl) = symba_plA%helio%swiftest%xh(:,1:npl) / rscale
         symba_plwksp%helio%swiftest%vh(:,1:npl) = symba_plA%helio%swiftest%vh(:,1:npl) / vscale
         symba_plwksp%helio%swiftest%xb(:,1:npl) = symba_plA%helio%swiftest%xb(:,1:npl) / rscale
         symba_plwksp%helio%swiftest%vb(:,1:npl) = symba_plA%helio%swiftest%vb(:,1:npl) / vscale
         symba_plwksp%helio%swiftest%rot(:,1:npl) = symba_plA%helio%swiftest%rot(:,1:npl) * tscale
         symba_plwksp%helio%swiftest%Ip(:,1:npl) = symba_plA%helio%swiftest%Ip(:,1:npl)
   
         if (linclude_fragments) then ! Append the fragments if they are included
            ! Energy calculation requires the fragments to be in the system barcyentric frame, s
            symba_plwksp%helio%swiftest%Ip(:,npl+1:npl_new) = Ip_frag(:,1:nfrag)
            symba_plwksp%helio%swiftest%mass(npl+1:npl_new) = m_frag(1:nfrag)
            symba_plwksp%helio%swiftest%radius(npl+1:npl_new) = rad_frag(1:nfrag)
            symba_plwksp%helio%swiftest%xb(:,npl+1:npl_new) =  xb_frag(:,1:nfrag)
            symba_plwksp%helio%swiftest%vb(:,npl+1:npl_new) =  vb_frag(:,1:nfrag)
            symba_plwksp%helio%swiftest%rot(:,npl+1:npl_new) = rot_frag(:,1:nfrag)
            symba_plwksp%helio%swiftest%status(npl+1:npl_new) = COLLISION
            call coord_b2h(npl_new, symba_plwksp%helio%swiftest)
            allocate(ltmp(npl_new))
            ltmp(1:npl) = lexclude(1:npl)
            ltmp(npl+1:npl_new) = .false.
            call move_alloc(ltmp, lexclude)
         end if
   
         where (lexclude(1:npl))
            symba_plwksp%helio%swiftest%status(1:npl) = INACTIVE
         end where
   
         nplm = count(symba_plwksp%helio%swiftest%mass > param%mtiny / mscale)
         call util_dist_index_plpl(npl_new, nplm, symba_plwksp)
         call symba_energy_eucl(npl_new, symba_plwksp, param%j2rp2, param%j4rp4, ke_orbit, ke_spin, pe, Lorbit, Lspin)
   
         ! Restore the big array
         deallocate(symba_plwksp%helio%swiftest%k_plpl) 
         nplm = count(symba_plA%helio%swiftest%mass > param%mtiny)
         if (lk_plpl) call util_dist_index_plpl(npl, nplm, symba_plA)

         ! Calculate the current fragment energy and momentum balances
         if (linclude_fragments) then
            Ltot_after(:) = Lorbit(:) + Lspin(:)
            Lmag_after = norm2(Ltot_after(:))
            ke_orbit_after = ke_orbit
            ke_spin_after = ke_spin
            pe_after = pe
            Etot_after = ke_orbit_after + ke_spin_after + pe_after
            dEtot = (Etot_after - Etot_before) / abs(Etot_before)
            dLmag = norm2(Ltot_after(:) - Ltot_before(:)) / Lmag_before
            ke_frag = 0._DP
            do i = 1, nfrag
               ke_frag = ke_frag + 0.5_DP * m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
            end do
            ke_target = ke_family + (ke_spin_before - ke_spin) + (pe_before - pe) - Qloss
            ke_offset = ke_frag - ke_target
            L_offset(:) = Ltot_before(:) - Ltot_after(:)
         else
            Ltot_before(:) = Lorbit(:) + Lspin(:)
            Lmag_before = norm2(Ltot_before(:))
            ke_orbit_before = ke_orbit
            ke_spin_before = ke_spin
            pe_before = pe
            Etot_before = ke_orbit_before + ke_spin_before + pe_before
         end if
      end associate
   return
   end subroutine calculate_system_energy

   subroutine shift_vector_to_origin(m_frag, vec_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the position or velocity of the fragments as needed to align them with the origin
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame

      ! Internals
      real(DP), dimension(NDIM)               :: mvec_frag, COM_offset
      integer(I4B)                            :: i

      mvec_frag(:) = 0.0_DP

      do i = 1, nfrag
         mvec_frag = mvec_frag(:) + vec_frag(:,i) * m_frag(i)
      end do
      COM_offset(:) = -mvec_frag(:) / mtot
      do i = 1, nfrag 
         vec_frag(:, i) = vec_frag(:, i) + COM_offset(:)
      end do

      return
   end subroutine shift_vector_to_origin

   subroutine set_fragment_position_vectors()
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
      !! pre-impact angular momentum. They are distributed on an ellipse surrounding the center of mass.
      !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.

      implicit none
      real(DP)   :: r_max, dis
      real(DP), dimension(NDIM) ::  L, L_sigma
      logical, dimension(:), allocatable :: loverlap
      integer(I4B) :: i, j


      allocate(loverlap(nfrag))

      ! Place the fragments into a region that is big enough that we should usually not have overlapping bodies
      ! An overlapping bodies will collide in the next time step, so it's not a major problem if they do (it just slows the run down)
      r_max = r_max_start

      ! We will treat the first fragment of the list as a special case. It gets initialized the maximum distance along the original impactor distance vector.
      ! This is done because in a regular disruption, the first body is the largest fragment.
      x_frag(:, 1) = -y_col_unit(:) * r_max  
      x_frag(:, 2) =  y_col_unit(:) * r_max  

      call random_number(x_frag(:,2:nfrag))
      loverlap(:) = .true.
      do while (any(loverlap(3:nfrag)))
         do i = 3, nfrag
            if (loverlap(i)) then
               call random_number(x_frag(:,i))
               x_frag(:, i) = 2 * (x_frag(:, i) - 0.5_DP) * r_max
            end if
         end do
         loverlap(:) = .false.
         do j = 1, nfrag
            do i = j + 1, nfrag
               dis = norm2(x_frag(:,j) - x_frag(:,i))
               loverlap(i) = loverlap(i) .or. (dis <= (rad_frag(i) + rad_frag(j))) 
            end do
         end do
         r_max = r_max + 0.2_DP
      end do
      call shift_vector_to_origin(m_frag, x_frag)

      do i = 1, nfrag
         rmag(i) = norm2(x_frag(:, i))
         v_r_unit(:, i) = x_frag(:, i) / rmag(i)
         call random_number(L_sigma(:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                        ! otherwise we can get an ill-conditioned system
         L(:) = z_col_unit(:) + 2e-3_DP * (L_sigma(:) - 0.5_DP)
         L(:) = L(:) / norm2(L(:)) 
         call util_crossproduct(L(:), v_r_unit(:, i), v_t_unit(:, i))
         xb_frag(:,i) = x_frag(:,i) + xcom(:)
      end do

      return
   end subroutine set_fragment_position_vectors

   subroutine set_fragment_tangential_velocities()
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
      implicit none
      ! Internals
      integer(I4B)                            :: i
      real(DP)                                :: L_orb_mag
      real(DP), dimension(2 * NDIM, 2 * NDIM) :: A ! LHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(2 * NDIM)           :: b  ! RHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(NDIM)               :: L_lin_others, L_orb_others, L

      v_frag(:,:) = 0.0_DP
      
      ! Divide up the pre-impact spin angular momentum equally between the various bodies by mass
      do i = 1, nfrag
         rot_frag(:,i) = L_frag_spin(:) / nfrag / (Ip_frag(:, i) * m_frag(i) * rad_frag(i)**2)
      end do

      L_orb_mag = norm2(L_frag_orb(:)) 
      ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
      ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
      ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
      L_lin_others(:) = 0.0_DP
      L_orb_others(:) = 0.0_DP
      do i = 1, nfrag
         if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
            A(1:3, i) = m_frag(i) * v_t_unit(:, i) 
            call util_crossproduct(v_r_unit(:, i), v_t_unit(:, i), L(:))
            A(4:6, i) = m_frag(i) * rmag(i) * L(:)
         else !For the remining bodies, distribute the angular momentum equally amongs them
            v_t_mag(i) = L_orb_mag / (m_frag(i) * rmag(i) * nfrag) 
            v_frag(:, i) = v_t_mag(i) * v_t_unit(:, i)
            L_lin_others(:) = L_lin_others(:) + m_frag(i) * v_frag(:, i)
            call util_crossproduct(x_frag(:, i), v_frag(:, i), L(:))
            L_orb_others(:) = L_orb_others(:) + m_frag(i) * L(:)
         end if
      end do
      b(1:3) = -L_lin_others(:)
      b(4:6) = L_frag_orb(:) - L_orb_others(:)
      v_t_mag(1:6) = util_solve_linear_system(A, b, 6, lmerge)
      if (lmerge) return
      do i = 1, 6
         v_frag(:, i) = v_t_mag(i) * v_t_unit(:, i)
      end do

      do i = 1, nfrag
         vb_frag(:,i) = v_frag(:,i) + vcom(:)
      end do

      return 
   end subroutine set_fragment_tangential_velocities

   subroutine set_fragment_radial_velocities(lmerge)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy.
      !! It will check that we don't end up with negative energy (bound system). If so, we'll set the fragment velocities to
      !! zero in the center of mass frame and indicate the the fragmentation should instead by a merger.
      !! It takes in the initial "guess" of velocities and solve for the a scaling factor applied to the radial component wrt the
      !! center of mass frame needed to correct the kinetic energy of the fragments in the system barycenter frame to match that of 
      !! the target kinetic energy required to satisfy the constraints.
      use symba_frag_pos_lambda_implementation
      implicit none
      ! Arguments
      logical,                intent(out)   :: lmerge
      ! Internals
      real(DP), parameter                   :: TOL = 1e-9_DP
      real(DP), dimension(:), allocatable   :: vflat 
      logical                               :: lerr
      integer(I4B)                          :: i
      real(DP), dimension(:), allocatable   :: v_r_initial, v_r_sigma
      real(DP), dimension(:,:), allocatable :: v_r

      ! Set the "target" ke_orbit_after (the value of the orbital kinetic energy that the fragments ought to have)
      
      if ((ke_target < 0.0_DP) .or. (ke_offset > 0.0_DP)) then
         !if (ke_target < 0.0_DP) write(*,*) 'Negative ke_target: ',ke_target
         !if (ke_offset > 0.0_DP) write(*,*) 'Positive ke_offset: ',ke_offset
         lmerge = .true.
         return
      end if

      allocate(v_r_initial, source=v_r_mag)
      ! Initialize radial velocity magnitudes with a random value that is approximately 10% of that found by distributing the kinetic energy equally
      allocate(v_r_sigma, source=v_r_mag)
      call random_number(v_r_sigma(1:nfrag))
      v_r_sigma(1:nfrag) = sqrt((1.0_DP + 2 * (v_r_sigma(1:nfrag) - 0.5_DP) * 1e-3_DP)) 
      v_r_initial(1:nfrag) = v_r_sigma(1:nfrag) * sqrt(abs(ke_offset) / (2 * m_frag(1:nfrag) * nfrag)) 

      ! Initialize the lambda function using a structure constructor that calls the init method
      ! Minimize the ke objective function using the BFGS optimizer
      v_r_mag(1:nfrag) = util_minimize_bfgs(ke_constraint(ke_objective_function, v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:), ke_target), &
                                             nfrag, v_r_initial(1:nfrag), TOL, lerr)
      if (lerr) then
         ! No solution exists for this case, so we need to indicate that this should be a merge
         ! This may happen due to setting the tangential velocities too high when setting the angular momentum constraint
         lmerge = .true.
         v_frag(:,:) = 0.0_DP
         do i = 1, nfrag
            vb_frag(:, i) = vcom(:)
         end do
      else
         lmerge = .false.
         ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
         vb_frag(:,1:nfrag) =  vmag_to_vb(v_r_mag(1:nfrag), v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:)) 
         do i = 1, nfrag
            v_frag(:, i) = vb_frag(:, i) - vcom(:)
         end do
      end if

      return
   end subroutine set_fragment_radial_velocities

   function ke_objective_function(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom, ke_target) result(fval) 
      !! Author: David A. Minton
      !!
      !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
      real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
      real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
      real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
      real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
      real(DP),                 intent(in)  :: ke_target !! Target kinetic energy
      ! Result
      real(DP)                              :: fval      !! The objective function result, which is the square of the difference between the calculated fragment kinetic energy and our target
                                                         !! Minimizing this brings us closer to our objective
      ! Internals
      integer(I4B)                        :: i
      real(DP), dimension(NDIM)           :: L
      real(DP), dimension(:,:), allocatable :: v_shift

      allocate(v_shift, mold=vb_frag)
      v_shift(:,:) = vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) 
      fval = -ke_target
      do i = 1, nfrag
         fval = fval + 0.5_DP * m_frag(i) * dot_product(v_shift(:, i), v_shift(:, i))
      end do
      ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
      fval = fval**2 

      return

   end function ke_objective_function

   function vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
      !! Author: David A. Minton
      !!
      !! Converts radial and tangential velocity magnitudes into barycentric velocity
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
      real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
      real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
      real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
      real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
      ! Result
      real(DP), dimension(:,:), allocatable   :: vb
      ! Internals
      integer(I4B) :: i

      allocate(vb, mold=v_r_unit)
      ! Make sure the velocity magnitude stays positive
      do i = 1, nfrag
         vb(:,i) = abs(v_r_mag(i)) * v_r_unit(:, i)
      end do
      ! In order to keep satisfying the kinetic energy constraint, we must shift the origin of the radial component of the velocities to the center of mass
      call shift_vector_to_origin(m_frag, vb)
      
      do i = 1, nfrag
         vb(:, i) = vb(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
      end do

   end function vmag_to_vb

   subroutine restructure_failed_fragments()
      !! Author: David A. Minton
      !!
      !! We failed to find a set of positions and velocities that satisfy all the constraints, and so we will alter the fragments and try again.
      !! Currently, this code will take the last fragment in the list and merge it with the first or second, alternating back and forth each time 
      !! it is called. 
      implicit none
      integer(I4B) :: i

      if (ke_offset < 0.0_DP) then 
         if (iflip == 1) then
            iflip = 2
         else
            iflip = 1
         end if
         m_frag(iflip) = m_frag(iflip) + m_frag(nfrag)
         rad_frag(iflip) = (rad_frag(iflip)**3 + rad_frag(nfrag)**3)**(1._DP / 3._DP) 
         m_frag(nfrag) = 0.0_DP
         rad_frag(nfrag) = 0.0_DP
         nfrag = nfrag - 1
      end if
      r_max_start = r_max_start + 1.0_DP ! The larger lever arm can help if the problem is in the angular momentum step
   end subroutine restructure_failed_fragments


end subroutine symba_frag_pos
