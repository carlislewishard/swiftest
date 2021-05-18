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
      real(DP), dimension(:,:), allocatable :: x_frag
      real(DP), dimension(NDIM)             :: L_target
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
      function abstract_objective_func(v_r_mag, v_r_unit, v_t_mag, v_t_unit, x_frag, m_frag, L_target, ke_target) result(fnorm)
         ! Template for the kinetic energy constraint function used for minimizing
         import DP
         real(DP), dimension(:),   intent(in) :: v_r_mag   !! Radial velocity mangitude
         real(DP), dimension(:,:), intent(in) :: v_r_unit  !! Radial velocity unit vector
         real(DP), dimension(:),   intent(in) :: v_t_mag   !! Tangential velocity magnitude
         real(DP), dimension(:,:), intent(in) :: v_t_unit  !! Tangential velocity unit vector
         real(DP), dimension(:,:), intent(in) :: x_frag    !! Position vectors
         real(DP), dimension(:),   intent(in) :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in) :: L_target  !! Target orbital momentum
         real(DP),                 intent(in) :: ke_target !! Target kinetic energ
         real(DP)                             :: fnorm     !! The objective function result: norm of the vector composed of the tangential momentum and energy
      end function
   end interface

   contains
      type(ke_constraint) function ke_objective_func_init(lambda, v_r_unit, v_t_mag, v_t_unit, x_frag, m_frag, L_target, ke_target)
         implicit none
         ! Arguments
         procedure(abstract_objective_func)   :: lambda    !! The lambda function
         real(DP), dimension(:,:), intent(in) :: v_r_unit  !! Radial velocity unit vector
         real(DP), dimension(:),   intent(in) :: v_t_mag   !! Tangential velocity magnitude
         real(DP), dimension(:,:), intent(in) :: v_t_unit  !! Tangential velocity unit vector
         real(DP), dimension(:,:), intent(in) :: x_frag    !! Position vectors
         real(DP), dimension(:),   intent(in) :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in) :: L_target  !! Target orbital momentum
         real(DP),                 intent(in) :: ke_target !! Target kinetic energ
         ! Internals
         associate(self => ke_objective_func_init)
            self%ke_objective_func_ptr  => lambda
            allocate(self%v_r_unit, source=v_r_unit)
            allocate(self%v_t_mag, source=v_t_mag)
            allocate(self%v_t_unit, source=v_t_unit)
            allocate(self%m_frag, source=m_frag)
            allocate(self%x_frag, source=x_frag)
            self%L_target(:) = L_target(:)
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
         if (allocated(self%x_frag)) deallocate(self%x_frag)
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
         fval = self%ke_objective_func_ptr(x, self%v_r_unit, self%v_t_mag, self%v_t_unit, self%x_frag, self%m_frag, self%L_target, self%ke_target)
      else
         error stop "KE Objective function was not initialized."
      end if
   end function ke_objective_func_eval

end module symba_frag_pos_lambda_implementation

subroutine symba_frag_pos (param, symba_plA, family, x, v, L_spin, Ip, mass, radius, &
                           Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, lmerge, Qloss)
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
   real(DP), intent(inout)                 :: Qloss
   real(DP), dimension(:,:), intent(inout) :: x, v, L_spin, Ip
   real(DP), dimension(:),   intent(inout) :: mass, radius, m_frag, rad_frag
   real(DP), dimension(:,:), intent(inout) :: Ip_frag
   real(DP), dimension(:,:), intent(inout) :: xb_frag, vb_frag, rot_frag
   logical, intent(out)                    :: lmerge ! Answers the question: Should this have been a merger instead?
   ! Internals
   real(DP), parameter                     :: f_spin = 0.20_DP !! Fraction of pre-impact orbital angular momentum that is converted to fragment spin
   real(DP)                                :: mscale, rscale, vscale, Lscale, tscale ! Scale factors that reduce quantities to O(~1) in the collisional system
   integer(I4B)                            :: i, j, nfrag, fam_size
   logical, dimension(:), allocatable      :: lexclude
   real(DP), dimension(NDIM, 2)            :: rot, L_orb 
   real(DP), dimension(:,:), allocatable   :: x_frag, v_frag, v_r_unit, v_t_unit
   real(DP), dimension(:), allocatable     :: rmag, v_r_mag, v_t_mag
   real(DP), dimension(NDIM)               :: xcom, vcom, Ltot_before, Ltot_after, L_residual
   real(DP), dimension(NDIM)               :: L_frag_spin, L_frag_tot, L_frag_orb
   real(DP)                                :: mtot, Lmag_before, Lmag_after
   real(DP)                                :: Etot_before, Etot_after, ke_orb_before, pe_before
   real(DP)                                :: pe_after, ke_spin_before, ke_spin_after, ke_orb_after, ke_family, ke_target, ke_frag
   real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
   character(len=*), parameter             :: fmtlabel = "(A14,10(ES9.2,1X,:))"

   allocate(x_frag, source=xb_frag)
   allocate(v_frag, source=vb_frag)
   nfrag = size(m_frag)

   associate(npl => symba_plA%helio%swiftest%nbody, status => symba_plA%helio%swiftest%status) 
      allocate(lexclude(npl))
      where (status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
         lexclude(1:npl) = .true.  
      elsewhere
         lexclude(1:npl) = .false. 
      end where
   end associate

   call set_scale_factors()
   call define_coordinate_system()

   call calculate_system_energy(ke_orb_before, ke_spin_before, pe_before, Ltot_before, linclude_fragments=.false.)
   Lmag_before = norm2(Ltot_before(:))
   Etot_before = ke_orb_before + ke_spin_before + pe_before

   call define_pre_collisional_family()
   call set_fragment_position_vectors()

   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,        "('              Energy normalized by |Etot_before|')")
   write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,        "('  First pass to get angular momentum ')")
   write(*,        "(' ---------------------------------------------------------------------------')")

   call set_fragment_tangential_velocities()
   
   call calculate_system_energy(ke_orb_after, ke_spin_after, pe_after, Ltot_after, linclude_fragments=.true.)
   Etot_after = ke_orb_after + ke_spin_after + pe_after
   Lmag_after = norm2(Ltot_after(:))

   write(*,fmtlabel) ' change      |',(ke_orb_after - ke_orb_before) / abs(Etot_before), &
                                       (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                       (ke_orb_after + ke_spin_after - ke_orb_before - ke_spin_before)/ abs(Etot_before), &
                                       (pe_after - pe_before) / abs(Etot_before), &
                                       (Etot_after - Etot_before) / abs(Etot_before), &
                                       norm2(Ltot_after - Ltot_before) / Lmag_before

   call set_fragment_radial_velocities(lmerge)
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,        "('  Second pass to get energy ')")
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,fmtlabel) ' Qloss      |',-Qloss / abs(Etot_before)
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,fmtlabel) ' T_family    |',ke_family / abs(Etot_before)
   write(*,fmtlabel) ' T_frag targ |',ke_target / abs(Etot_before)
   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,fmtlabel) ' T_frag calc |',ke_frag / abs(Etot_before)
   write(*,fmtlabel) ' residual    |',(ke_frag - ke_target) / abs(Etot_before)

   call calculate_system_energy(ke_orb_after, ke_spin_after, pe_after, Ltot_after, linclude_fragments=.true.)
   Etot_after = ke_orb_after + ke_spin_after + pe_after
   Lmag_after = norm2(Ltot_after(:))

   write(*,        "(' ---------------------------------------------------------------------------')")
   write(*,fmtlabel) ' change      |',(ke_orb_after - ke_orb_before) / abs(Etot_before), &
                                       (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                       (ke_orb_after + ke_spin_after - ke_orb_before - ke_spin_before)/ abs(Etot_before), &
                                       (pe_after - pe_before) / abs(Etot_before), &
                                       (Etot_after - Etot_before) / abs(Etot_before), &
                                       norm2(Ltot_after - Ltot_before) / Lmag_before

   lmerge = lmerge .or. ((Etot_after - Etot_before) / abs(Etot_before) > 0._DP) 

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

      mtot = 1.0_DP

      ! Set scale factors
      mscale = sum(mass(:)) 
      rscale = norm2(x(:,2) - x(:,1))
      tscale = rscale / norm2(v(:,2) - v(:,1))
      vscale = rscale / tscale
      Lscale = mscale * rscale * vscale

      mass = mass / mscale
      radius = radius / rscale
      x = x / rscale
      v = v / vscale
      L_spin = L_spin / Lscale

      m_frag = m_frag / mscale
      rad_frag = rad_frag / rscale
      Qloss = Qloss / (mscale * vscale**2) 
      return
   end subroutine set_scale_factors

   subroutine restore_scale_factors()
      !! author: David A. Minton
      !!
      !! Restores dimenional quantities back to the system units
      implicit none
      mass = mass * mscale
      radius = radius * rscale
      x = x * rscale
      v = v * vscale
      L_spin = L_spin * Lscale

      xb_frag = xb_frag * rscale
      vb_frag = vb_frag * vscale
      m_frag = m_frag * mscale
      rot_frag = rot_frag / tscale
      rad_frag = rad_frag * rscale
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

      ! Find the center of mass of the collisional system	
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

      L_orb(:, :) = 0.0_DP
      ! Compute orbital angular momentum of pre-impact system
      do j = 1, 2
         xc(:) = x(:, j) - xcom(:) 
         vc(:) = v(:, j) - vcom(:)
         call util_crossproduct(xc(:), vc(:), x_cross_v(:))
         L_orb(:, j) = mass(j) * x_cross_v(:)
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
      associate(vbpl => symba_plA%helio%swiftest%vb, Mpl => symba_plA%helio%swiftest%mass)
         fam_size = size(family)

         ! We need the original kinetic energy of just the pre-impact family members in order to balance the energy later
         ke_family = 0.0_DP
         do i = 1, fam_size
            ke_family = ke_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i)))
            lexclude(family(i)) = .true. ! For all subsequent energy calculations the pre-impact family members will be replaced by the fragments
         end do
         ke_family = 0.5_DP * ke_family / (mscale * vscale**2)
      end associate
      return
   end subroutine define_pre_collisional_family

   subroutine calculate_system_energy(ke_orb, ke_spin, pe, Ltot, linclude_fragments)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy, including all bodies in the symba_plA list that do not have a corresponding value of the lexclude array that is true
      !! and optionally including fragments.
      use module_swiftestalloc
      implicit none
      ! Arguments
      real(DP),               intent(out) :: ke_orb, ke_spin, pe
      real(DP), dimension(:), intent(out) :: Ltot
      logical,                intent(in) :: linclude_fragments
      ! Internals
      integer(I4B) :: i, npl_new, nplm
      logical, dimension(:), allocatable :: ltmp
      logical :: lk_plpl
      real(DP) :: te
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
         symba_plwksp%helio%swiftest%rhill(1:npl) = symba_plA%helio%swiftest%rhill(1:npl) / rscale
         symba_plwksp%helio%swiftest%xb(:,1:npl) = symba_plA%helio%swiftest%xb(:,1:npl) / rscale
         symba_plwksp%helio%swiftest%vb(:,1:npl) = symba_plA%helio%swiftest%vb(:,1:npl) / vscale
         symba_plwksp%helio%swiftest%rot(:,1:npl) = symba_plA%helio%swiftest%rot(:,1:npl) * tscale
         symba_plwksp%helio%swiftest%Ip(:,1:npl) = symba_plA%helio%swiftest%Ip(:,1:npl)
   
         if (linclude_fragments) then ! Append the fragments if they are included
            ! Energy calculation requires the fragments to be in the system barcyentric frame, s
            symba_plwksp%helio%swiftest%Ip(:,npl+1:npl_new) = Ip_frag(:,:)
            symba_plwksp%helio%swiftest%mass(npl+1:npl_new) = m_frag(:)
            symba_plwksp%helio%swiftest%radius(npl+1:npl_new) = rad_frag(:)
            symba_plwksp%helio%swiftest%xb(:,npl+1:npl_new) =  xb_frag(:,:)
            symba_plwksp%helio%swiftest%vb(:,npl+1:npl_new) =  vb_frag(:,:)
            symba_plwksp%helio%swiftest%rot(:,npl+1:npl_new) = rot_frag(:,:)
            symba_plwksp%helio%swiftest%status(npl+1:npl_new) = COLLISION
            call coord_b2h(npl_new, symba_plwksp%helio%swiftest)
            allocate(ltmp(npl_new))
            ltmp(1:npl) = lexclude(1:npl)
            ltmp(npl+1:npl_new) = .false.
            call move_alloc(ltmp, lexclude)
   
            ke_frag = 0._DP
            do i = 1, nfrag
               ke_frag = ke_frag + 0.5_DP * m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
            end do
         end if
   
         where (lexclude(1:npl))
            symba_plwksp%helio%swiftest%status(1:npl) = INACTIVE
         end where
   
         nplm = count(symba_plwksp%helio%swiftest%mass > param%mtiny)
         call util_dist_index_plpl(npl_new, nplm, symba_plwksp)
         call symba_energy_eucl(npl_new, symba_plwksp, param%j2rp2, param%j4rp4, ke_orb, ke_spin, pe, te, Ltot)
   
         ! Restore the big array
         deallocate(symba_plwksp%helio%swiftest%k_plpl) 
         nplm = count(symba_plA%helio%swiftest%mass > param%mtiny)
         if (lk_plpl) call util_dist_index_plpl(npl, nplm, symba_plA)
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
      integer(I4B)                            :: i, nfrag
      real(DP)                                :: mtot

      nfrag = size(m_frag)
      mtot = sum(m_frag)
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
      real(DP)   :: r_max
      real(DP), dimension(NDIM) ::  L, L_sigma

      allocate(rmag(nfrag))
      allocate(v_r_mag(nfrag))
      allocate(v_t_mag(nfrag))
      allocate(v_r_unit(NDIM,nfrag))
      allocate(v_t_unit(NDIM,nfrag))

      ! Re-normalize position and velocity vectors by the fragment number so that for our initial guess we weight each
      ! fragment position by the mass and assume equipartition of energy for the velocity
      r_max = 1.5_DP * sum(rad_frag(:)) / PI

      ! We will treat the first fragment of the list as a special case.
      x_frag(:, 1) = -z_col_unit(:) 
      call random_number(x_frag(:,2:nfrag)) 
      
      x_frag(:, :) = x_frag(:, :) * r_max
      call shift_vector_to_origin(m_frag, x_frag)

      do i = 1, nfrag
         xb_frag(:,i) = x_frag(:,i) + xcom(:)
         rmag(i) = norm2(x_frag(:, i))
         v_r_unit(:, i) = x_frag(:, i) / rmag(i)
         call random_number(L_sigma(:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                        ! otherwise we can get an ill-conditioned system
         L(:) = z_col_unit(:) + 2e-3_DP * (L_sigma(:) - 0.5_DP)
         L(:) = L(:) / norm2(L(:)) 
         call util_crossproduct(L(:), v_r_unit(:, i), v_t_unit(:, i))
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
      logical,                  intent(out)   :: lmerge
      ! Internals
      real(DP), parameter                   :: TOL = epsilon(1._DP)
      real(DP), dimension(:), allocatable   :: vflat 
      logical                               :: lerr
      integer(I4B)                          :: i
      real(DP)                              :: ke_tangential
      real(DP), dimension(:), allocatable   :: v_r_initial
      real(DP), dimension(:,:), allocatable :: v_r

      ! Set the "target" ke_orb_after (the value of the orbital kinetic energy that the fragments ought to have)
      ke_target = ke_family + (ke_spin_before - ke_spin_after) + (pe_before - pe_after) - Qloss
      if (ke_target < 0.0_DP) then
         lmerge = .true.
         return
      end if

      ! Initialize radial velocity magnitudes with a random value:
      allocate(v_r_initial, source=v_r_mag)
      call random_number(v_r_initial(:))
      v_r_initial(:) = 3 * v_r_initial(:) * sqrt(2 * ke_target / nfrag / m_frag(:)) 
      ! Initialize the lambda function using a structure constructor that calls the init method
      ! Minimize error using the BFGS optimizer
      v_r_mag(:) = util_minimize_bfgs(ke_constraint(ke_objective_function, v_r_unit, v_t_mag, v_t_unit, x_frag, m_frag, L_frag_orb, ke_target), nfrag, v_r_initial, TOL, lerr)
      if (lerr) then
         ! No solution exists for this case, so we need to indicate that this should be a merge
         ! This may happen due to setting the tangential velocities too high when setting the angular momentum constraint
         lmerge = .true.
         v_frag(:,:) = 0.0_DP
      else
         lmerge = .false.
         ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
         allocate(v_r, mold=v_frag)
         do i = 1, nfrag
            v_r(:, i) = v_r_mag(i) * v_r_unit(:, i)
         end do
         call shift_vector_to_origin(m_frag, v_r)
        
         ! Recombine the tangential and radial components into the final velocity vector
         do i = 1, nfrag
            v_r_mag(i) = dot_product(v_r(:,i), v_r_unit(:, i))
            v_frag(:, i) = v_r_mag(i) * v_r_unit(:, i) + v_t_mag(i) * v_t_unit(:, i)
         end do
         call shift_vector_to_origin(m_frag, v_frag)
      end if

      do i = 1, nfrag
         vb_frag(:,i) = v_frag(:,i) + vcom(:)
      end do

      return
   end subroutine set_fragment_radial_velocities

   function ke_objective_function(v_r_mag, v_r_unit, v_t_mag, v_t_unit, x_frag, m_frag, L_target, ke_target) result(fval) 
      ! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
      real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
      real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
      real(DP), dimension(:,:), intent(in)  :: x_frag    !! Velocity and position vectors
      real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
      real(DP), dimension(:),   intent(in)  :: L_target  !! Target orbital momentum
      real(DP),                 intent(in)  :: ke_target !! Target kinetic energ
      ! Result
      real(DP)                              :: fval           !! The objective function result: norm of the vector composed of the tangential momentum and energy
                                                                          !! Minimizing this brings us closer to our objective
      ! Internals
      integer(I4B)                        :: i, nfrag, nsol
      real(DP), dimension(NDIM)           :: L
      real(DP), dimension(:,:), allocatable :: v_shift

      nfrag = size(m_frag)
      allocate(v_shift, mold=v_r_unit)
      ! In order to keep satisfying the kinetic energy constraint, we must shift the origin of the radial component of the velocities to the center of mass
      do i = 1, nfrag
         v_shift(:,i) = v_r_mag(i) * v_r_unit(:, i)
      end do
      call shift_vector_to_origin(m_frag, v_shift)
      
      fval = -ke_target
      do i = 1, nfrag
         v_shift(:, i) = v_shift(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
         fval = fval + 0.5_DP * m_frag(i) * dot_product(v_shift(:, i), v_shift(:, i))
      end do
      fval = fval**4

      return

   end function ke_objective_function


end subroutine symba_frag_pos
