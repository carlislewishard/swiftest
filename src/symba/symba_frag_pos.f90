subroutine symba_frag_pos(param, symba_plA, family, x, v, L_spin, Ip, mass, radius, &
                          nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lfailure)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Places the collision fragments on a circle oriented with a plane defined
   !! by the position and velocity vectors of the collision
   !! 
   use, intrinsic :: ieee_exceptions
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use lambda_function
   use module_interfaces, except_this_one => symba_frag_pos
   implicit none
   ! Arguments
   type(swiftest_parameters), intent(in) :: param 
   type(symba_pl), intent(inout)           :: symba_plA
   integer(I4B), dimension(:), intent(in)  :: family
   real(DP), dimension(:,:), intent(inout) :: x, v, L_spin, Ip
   real(DP), dimension(:),   intent(inout) :: mass, radius
   integer(I4B), intent(inout)             :: nfrag
   real(DP), dimension(:), allocatable,   intent(inout) :: m_frag, rad_frag
   real(DP), dimension(:,:), allocatable, intent(inout) :: Ip_frag
   real(DP), dimension(:,:), allocatable, intent(inout) :: xb_frag, vb_frag, rot_frag
   logical, intent(out)                    :: lfailure ! Answers the question: Should this have been a merger instead?
   real(DP), intent(inout)                 :: Qloss
   ! Internals
   real(DP)                                :: mscale, rscale, vscale, tscale, Lscale, Escale ! Scale factors that reduce quantities to O(~1) in the collisional system
   real(DP)                                :: mtot 
   real(DP), dimension(NDIM)               :: xcom, vcom
   integer(I4B)                            :: ii
   logical, dimension(:), allocatable      :: lexclude
   real(DP), dimension(NDIM, 2)            :: rot, L_orb 
   real(DP), dimension(:,:), allocatable   :: x_frag, v_frag, v_r_unit, v_t_unit, v_h_unit
   real(DP), dimension(:), allocatable     :: rmag, rotmag, v_r_mag, v_t_mag
   real(DP), dimension(NDIM)               :: Ltot_before
   real(DP), dimension(NDIM)               :: Ltot_after
   real(DP)                                :: Etot_before, ke_orbit_before, ke_spin_before, pe_before, Lmag_before
   real(DP)                                :: Etot_after,  ke_orbit_after,  ke_spin_after,  pe_after,  Lmag_after, dEtot, dLmag
   real(DP), dimension(NDIM)               :: L_frag_tot, L_frag_orb
   real(DP)                                :: ke_frag_budget, ke_frag_orbit, ke_radial, ke_frag_spin, ke_avg_deficit, ke_avg_deficit_old
   real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
   character(len=*), parameter             :: fmtlabel = "(A14,10(ES11.4,1X,:))"
   integer(I4B)                            :: try, subtry
   integer(I4B), parameter                 :: NFRAG_MIN = 7 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
   real(DP)                                :: r_max_start, r_max_start_old, r_max, f_spin 
   real(DP), parameter                     :: Ltol = 10 * epsilon(1.0_DP)
   real(DP), parameter                     :: Etol = 1e-10_DP
   integer(I4B), parameter                 :: MAXTRY = 3000
   integer(I4B), parameter                 :: TANTRY = 3
   logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes, fpe_quiet_modes

   if (nfrag < NFRAG_MIN) then
      write(*,*) "symba_frag_pos needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
      lfailure = .true.
      return
   end if

   call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
   fpe_quiet_modes(:) = .false.
   call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

   f_spin = 0.05_DP
   mscale = 1.0_DP
   rscale = 1.0_DP
   vscale = 1.0_DP
   tscale = 1.0_DP
   Lscale = 1.0_DP
   Escale = 1.0_DP

   associate(npl => symba_plA%helio%swiftest%nbody, status => symba_plA%helio%swiftest%status) 
      allocate(lexclude(npl))
      where (status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
         lexclude(1:npl) = .true.  
      elsewhere
         lexclude(1:npl) = .false. 
      end where
   end associate

   allocate(x_frag, source=xb_frag)
   allocate(v_frag, source=vb_frag)

   call set_scale_factors()
   call define_coordinate_system()
   call calculate_system_energy(linclude_fragments=.false.)
  
   r_max_start = norm2(x(:,2) - x(:,1))
   try = 1
   lfailure = .false.
   ke_avg_deficit = 0.0_DP
   do while (try < MAXTRY)
      lfailure = .false.
      ke_avg_deficit_old = ke_avg_deficit
      ke_avg_deficit = 0.0_DP
      subtry = 1
      do 
         call set_fragment_position_vectors()
         call set_fragment_tangential_velocities(lfailure)
         ke_avg_deficit = ke_avg_deficit - ke_radial
         subtry = subtry + 1
         if (.not.lfailure .or. subtry == TANTRY) exit
         !write(*,*) 'Trying new arrangement'
      end do
      ke_avg_deficit = ke_avg_deficit / subtry
      if (lfailure) write(*,*) 'Failed to find tangential velocities'

      if (.not.lfailure) then
         call set_fragment_radial_velocities(lfailure)
         if (lfailure) write(*,*) 'Failed to find radial velocities'
         if (.not.lfailure) then
            call calculate_system_energy(linclude_fragments=.true.)
            !write(*,*) 'Qloss : ',Qloss
            !write(*,*) '-dEtot: ',-dEtot
            !write(*,*) 'delta : ',abs((dEtot + Qloss)) 
            if ((abs(dEtot + Qloss) > Etol) .or. (dEtot > 0.0_DP)) then
               write(*,*) 'Failed due to high energy error: ',dEtot, abs(dEtot + Qloss) / Etol
               lfailure = .true.
            else if (abs(dLmag) / Lmag_before > Ltol) then
               write(*,*) 'Failed due to high angular momentum error: ', dLmag / Lmag_before
               lfailure = .true.
            end if
         end if
      end if
 
      if (.not.lfailure) exit
      call restructure_failed_fragments()
      try = try + 1
   end do
   write(*,        "(' -------------------------------------------------------------------------------------')")
   write(*,        "('  Final diagnostic')")
   write(*,        "(' -------------------------------------------------------------------------------------')")
   if (lfailure) then
      write(*,*) "symba_frag_pos failed after: ",try," tries"
      do ii = 1, nfrag
         vb_frag(:, ii) = vcom(:)
      end do
   else
      write(*,*) "symba_frag_pos succeeded after: ",try," tries"
      write(*,        "(' dL_tot should be very small' )")
      write(*,fmtlabel) ' dL_tot      |', dLmag / Lmag_before
      write(*,        "(' dE_tot should be negative and equal to Qloss' )")
      write(*,fmtlabel) ' dE_tot      |', dEtot / abs(Etot_before)
      write(*,fmtlabel) ' Qloss       |', -Qloss / abs(Etot_before)
      write(*,fmtlabel) ' dE - Qloss  |', (Etot_after - Etot_before + Qloss) / abs(Etot_before)
   end if
   write(*,        "(' -------------------------------------------------------------------------------------')")

   call restore_scale_factors()
   call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily

   return 

   contains

   ! Because of the complexity of this procedure, we have chosen to break it up into a series of nested subroutines.

   subroutine set_scale_factors()
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. This scaling makes it easier for the non-linear minimization 
      !! to converge on a solution
      implicit none
      integer(I4B) :: i

      ! Find the center of mass of the collisional system	
      mtot = sum(mass(:)) 
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

      ! Set scale factors
      !! Because of the implied G, mass is actually G*mass with units of distance**3 / time**2
      Escale = 0.5_DP * (mass(1) * dot_product(v(:,1), v(:,1)) + mass(2) * dot_product(v(:,2), v(:,2)))
      rscale = sum(radius(:))
      mscale = sqrt(Escale * rscale) 
      vscale = sqrt(Escale / mscale) 
      tscale = rscale / vscale 
      Lscale = mscale * rscale * vscale

      xcom(:) = xcom(:) / rscale
      vcom(:) = vcom(:) / vscale

      mtot = mtot / mscale
      mass = mass / mscale
      radius = radius / rscale
      x = x / rscale
      v = v / vscale
      L_spin = L_spin / Lscale
      do i = 1, 2
         rot(:,i) = L_spin(:,i) / (mass(i) * radius(i)**2 * Ip(3, i))
      end do

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

      call ieee_set_halting_mode(IEEE_ALL,.false.)
      ! Restore scale factors
      xcom(:) = xcom(:) * rscale
      vcom(:) = vcom(:) * vscale

      mtot = mtot * mscale
      mass = mass * mscale
      radius = radius * rscale
      x = x * rscale
      v = v * vscale
      L_spin = L_spin * Lscale
      do i = 1, 2
         rot(:,i) = L_spin(:,i) * (mass(i) * radius(i)**2 * Ip(3, i))
      end do

      m_frag = m_frag * mscale
      rad_frag = rad_frag * rscale
      rot_frag = rot_frag / tscale
      x_frag = x_frag * rscale
      v_frag = v_frag * vscale
      Qloss = Qloss * Escale

      do i = 1, nfrag
         xb_frag(:, i) = x_frag(:, i) + xcom(:)
         vb_frag(:, i) = v_frag(:, i) + vcom(:)
      end do

      Etot_before = Etot_before * Escale
      pe_before = pe_before * Escale
      ke_spin_before = ke_spin_before * Escale
      ke_orbit_before = ke_orbit_before * Escale
      Ltot_before = Ltot_before * Lscale
      Lmag_before = Lmag_before * Lscale 
      Etot_after = Etot_after * Escale
      pe_after = pe_after * Escale
      ke_spin_after = ke_spin_after * Escale
      ke_orbit_after = ke_orbit_after * Escale
      Ltot_after = Ltot_after * Lscale
      Lmag_after = Lmag_after * Lscale 

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
      allocate(rotmag(nfrag))
      allocate(v_r_mag(nfrag))
      allocate(v_t_mag(nfrag))
      allocate(v_r_unit(NDIM,nfrag))
      allocate(v_t_unit(NDIM,nfrag))
      allocate(v_h_unit(NDIM,nfrag))

      rmag(:) = 0.0_DP
      rotmag(:) = 0.0_DP
      v_r_mag(:) = 0.0_DP
      v_t_mag(:) = 0.0_DP
      v_r_unit(:,:) = 0.0_DP
      v_t_unit(:,:) = 0.0_DP
      v_h_unit(:,:) = 0.0_DP

      L_orb(:, :) = 0.0_DP
      ! Compute orbital angular momentum of pre-impact system
      do i = 1, 2
         xc(:) = x(:, i) - xcom(:) 
         vc(:) = v(:, i) - vcom(:)
         call util_crossproduct(xc(:), vc(:), x_cross_v(:))
         L_orb(:, i) = mass(i) * x_cross_v(:)
      end do

      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      L_frag_tot(:) = L_spin(:, 1) + L_spin(:, 2) + L_orb(:, 1) + L_orb(:, 2)

      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))     
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))

      ! We will initialize fragments on a plane defined by the pre-impact system, with the z-axis aligned with the angular momentum vector
      ! and the y-axis aligned with the pre-impact distance vector.
      y_col_unit(:) = delta_r(:) / r_col_norm  
      z_col_unit(:) = L_frag_tot(:) / norm2(L_frag_tot)
      ! The cross product of the y- by z-axis will give us the x-axis
      call util_crossproduct(y_col_unit, z_col_unit, x_col_unit)

      return
   end subroutine define_coordinate_system

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
            lexclude(family(:)) = .true. 
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
   
         where (lexclude(1:npl_new))
            symba_plwksp%helio%swiftest%status(1:npl_new) = INACTIVE
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
            dEtot = Etot_after - Etot_before 
            dLmag = norm2(Ltot_after(:) - Ltot_before(:)) 
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
      real(DP)  :: dis, rad
      real(DP), dimension(NDIM) ::  L_sigma
      logical, dimension(:), allocatable :: loverlap
      integer(I4B) :: i, j

      allocate(loverlap(nfrag))

      ! Place the fragments into a region that is big enough that we should usually not have overlapping bodies
      ! An overlapping bodies will collide in the next time step, so it's not a major problem if they do (it just slows the run down)
      r_max = r_max_start
      rad = sum(radius(:))

      ! We will treat the first two fragments of the list as special cases. They get initialized the maximum distances apart along the original impactor distance vector.
      ! This is done because in a regular disruption, the first body is the largest, the second the second largest, and the rest are smaller equal-mass fragments.

      call random_number(x_frag(:,3:nfrag))
      loverlap(:) = .true.
      do while (any(loverlap(3:nfrag)))
         x_frag(:, 1) = x(:, 1) - xcom(:) 
         x_frag(:, 2) = x(:, 2) - xcom(:) 
         r_max = r_max + 0.1_DP * rad
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
      end do
      call shift_vector_to_origin(m_frag, x_frag)

      do i = 1, nfrag
         rmag(i) = norm2(x_frag(:, i))
         v_r_unit(:, i) = x_frag(:, i) / rmag(i)
         call random_number(L_sigma(:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                        ! otherwise we can get an ill-conditioned system
         v_h_unit(:, i) = z_col_unit(:) + 2e-1_DP * (L_sigma(:) - 0.5_DP)
         v_h_unit(:, i) = v_h_unit(:, i) / norm2(v_h_unit(:, i)) 
         call util_crossproduct(v_h_unit(:, i), v_r_unit(:, i), v_t_unit(:, i))
         xb_frag(:,i) = x_frag(:,i) + xcom(:)
      end do

      return
   end subroutine set_fragment_position_vectors

   subroutine set_fragment_tangential_velocities(lerr)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the tangential velocities and spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
      !! This procedure works in several stages, with a goal to solve the angular and linear momentum constraints on the fragments, while still leaving a positive balance of
      !! our fragment kinetic energy (ke_frag_budget) that we can put into the radial velocity distribution.
      !!
      !! The first thing we'll try to do is solve for the tangential velocities of the first 6 fragments, using angular and linear momentum as constraints and an initial
      !! tangential velocity distribution for the remaining bodies (if there are any) that distributes their angular momentum equally between them.
      !! If that doesn't work and we blow our kinetic energy budget, we will attempt to find a tangential velocity distribution that minimizes the kinetic energy while
      !! conserving momentum. 
      !!
      !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
      implicit none
      ! Arguments
      logical, intent(out)  :: lerr
      ! Internals
      integer(I4B) :: i
      real(DP), parameter                  :: TOL = 1e-4_DP
      real(DP), dimension(:), allocatable  :: v_t_initial
      type(lambda_obj)                     :: spinfunc
      type(lambda_obj_err)                 :: objective_function
      real(DP), dimension(NDIM) :: L_frag_spin, L_remainder, Li, rot_L, rot_ke

      ! Initialize the fragments with 0 velocity and spin so we can divide up the balance between the tangential, radial, and spin components while conserving momentum
      lerr = .false.
      vb_frag(:,:) = 0.0_DP
      rot_frag(:,:) = 0.0_DP
      v_t_mag(:) = 0.0_DP
      v_r_mag(:) = 0.0_DP

      call calculate_system_energy(linclude_fragments=.true.)
      ke_frag_budget = -dEtot - Qloss
      !write(*,*) '***************************************************'
      !write(*,*) 'Original dis   : ',norm2(x(:,2) - x(:,1))
      !write(*,*) 'r_max          : ',r_max
      !write(*,*) 'f_spin         : ',f_spin
      !write(*,*) '***************************************************'
      !write(*,*) 'Energy balance so far: '
      !write(*,*) 'ke_frag_budget : ',ke_frag_budget
      !write(*,*) 'ke_orbit_before: ',ke_orbit_before 
      !write(*,*) 'ke_orbit_after : ',ke_orbit_after  
      !write(*,*) 'ke_spin_before : ',ke_spin_before 
      !write(*,*) 'ke_spin_after  : ',ke_spin_after  
      !write(*,*) 'pe_before      : ',pe_before 
      !write(*,*) 'pe_after       : ',pe_after  
      !write(*,*) 'Qloss          : ',Qloss
      !write(*,*) '***************************************************'
      if (ke_frag_budget < 0.0_DP) then
         write(*,*) 'Negative ke_frag_budget: ',ke_frag_budget
         r_max_start = r_max_start / 2 
         lerr = .true.
         return
      end if

      allocate(v_t_initial, mold=v_t_mag)

      L_frag_spin(:) = 0.0_DP
      ke_frag_spin = 0.0_DP
      ! Start the first two bodies with the same rotation as the original two impactors, then distribute the remaining angular momentum among the rest
      do i = 1, 2
         rot_frag(:, i) = rot(:, i)
         L_frag_spin(:) = L_frag_spin(:) + m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i) * rot_frag(:, i)
      end do
      L_frag_orb(:) =  L_frag_tot(:) - L_frag_spin(:)
      L_frag_spin(:) = 0.0_DP
      do i = 1, nfrag
         ! Convert a fraction (f_spin) of either the remaining angular momentum or kinetic energy budget into spin, whichever gives the smaller rotation so as not to blow any budgets
         rot_ke(:) = sqrt(2 * f_spin * ke_frag_budget / (nfrag * m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i))) * L_frag_orb(:) / norm2(L_frag_orb(:))
         rot_L(:) = f_spin * L_frag_orb(:) / (nfrag * m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i))
         if (norm2(rot_ke) < norm2(rot_L)) then
            rot_frag(:,i) = rot_frag(:, i) + rot_ke(:)
         else
            rot_frag(:, i) = rot_frag(:, i) + rot_L(:)
         end if
         L_frag_spin(:) = L_frag_spin(:) + m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i) * rot_frag(:, i)
         ke_frag_spin = ke_frag_spin + m_frag(i) * Ip_frag(3, i) * rad_frag(i)**2 * dot_product(rot_frag(:, i), rot_frag(:, i))
      end do
      ke_frag_spin = 0.5_DP * ke_frag_spin
      ! Convert a fraction of the pre-impact angular momentum into fragment spin angular momentum
      L_frag_orb(:) =  L_frag_tot(:) - L_frag_spin(:)
      L_remainder(:) = L_frag_orb(:)
      ! Next we will solve for the tangential component of the velocities that both conserves linear momentum and uses the remaining angular momentum not used in spin.
      ! This will be done using a linear solver that solves for the tangential velocities of the first 6 fragments, constrained by the linear and angular momentum vectors, 
      ! which is embedded in a non-linear minimizer that will adjust the tangential velocities of the remaining i>6 fragments to minimize kinetic energy for a given momentum solution
      ! The initial conditions fed to the minimizer for the fragments will be the remaining angular momentum distributed between the fragments.
      do i = 1, nfrag
         v_t_initial(i) = norm2(L_remainder(:)) / ((nfrag - i + 1) * m_frag(i) * norm2(x_frag(:,i)))
         call util_crossproduct(m_frag(i) * x_frag(:,i), v_t_initial(i) * v_t_unit(:, i), Li(:))
         L_remainder(:) = L_remainder(:) - Li(:)
      end do

      ! Find the local kinetic energy minimum for the system that conserves linear and angular momentum
      objective_function = lambda_obj(tangential_objective_function, lerr)
      v_t_mag(7:nfrag) = util_minimize_bfgs(objective_function, nfrag-6, v_t_initial(7:nfrag), TOL, lerr)
      ! Now that the KE-minimized values of the i>6 fragments are found, calculate the momentum-conserving solution for tangential velociteis
      v_t_initial(7:nfrag) = v_t_mag(7:nfrag)
      v_t_mag(1:nfrag) = solve_fragment_tangential_velocities(v_t_mag_input=v_t_initial(7:nfrag), lerr=lerr)

      ! Perform one final shift of the radial velocity vectors to align with the center of mass of the collisional system (the origin)
      vb_frag(:,1:nfrag) = vmag_to_vb(v_r_mag(1:nfrag), v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:)) 
      ! Now do a kinetic energy budget check to make sure we are still within the budget.
      ke_frag_orbit = 0.0_DP
      do i = 1, nfrag
         v_frag(:, i) = vb_frag(:, i) - vcom(:)
         ke_frag_orbit = ke_frag_orbit + m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
      end do
      ke_frag_orbit = 0.5_DP * ke_frag_orbit
      ke_radial = ke_frag_budget - ke_frag_orbit - ke_frag_spin

      ! If we are over the energy budget, flag this as a failure so we can try again
      lerr = (ke_radial < 0.0_DP)
      !write(*,*) 'Tangential'
      !write(*,*) 'ke_frag_budget: ',ke_frag_budget
      !write(*,*) 'ke_frag_orbit : ',ke_frag_orbit
      !write(*,*) 'ke_frag_spin  : ',ke_frag_spin
      !write(*,*) 'ke_radial     : ',ke_radial

      return

   end subroutine set_fragment_tangential_velocities

   function tangential_objective_function(v_t_mag_input, lerr) result(fval)
      !! Author: David A. Minton
      !!
      !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_t_mag_input   !! Unknown tangential component of velocity vector set previously by angular momentum constraint
      logical,                  intent(out) :: lerr            !! Error flag
      ! Result
      real(DP)                              :: fval
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(:,:), allocatable :: v_shift
      real(DP), dimension(:), allocatable :: v_t_new
      real(DP) :: keo

      lerr = .false.

      allocate(v_shift(NDIM, nfrag))
      allocate(v_t_new(nfrag))

      v_t_new(:) = solve_fragment_tangential_velocities(v_t_mag_input=v_t_mag_input(:), lerr=lerr)
      v_shift(:,:) = vmag_to_vb(v_r_mag, v_r_unit, v_t_new, v_t_unit, m_frag, vcom) 

      keo = 0.0_DP
      do i = 1, nfrag
         keo = keo + m_frag(i) * dot_product(v_shift(:, i), v_shift(:, i))
      end do
      keo = 0.5_DP * keo
      fval = keo 
      lerr = .false.
      return
   end function tangential_objective_function

   function solve_fragment_tangential_velocities(lerr, v_t_mag_input) result(v_t_mag_output)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
      implicit none
      ! Arguments
      logical,                          intent(out) :: lerr            !! Error flag 
      real(DP), dimension(:), optional, intent(in)  :: v_t_mag_input   !! Unknown tangential velocities for fragments 7:nfrag
      ! Internals
      integer(I4B)                            :: i
      ! Result
      real(DP), dimension(:), allocatable     :: v_t_mag_output

      real(DP), dimension(2 * NDIM, 2 * NDIM) :: A ! LHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(2 * NDIM)           :: b  ! RHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(NDIM)               :: L_lin_others, L_orb_others, L, vtmp

      v_frag(:,:) = 0.0_DP
      lerr = .false.

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
         else if (present(v_t_mag_input)) then
            vtmp(:) = v_t_mag_input(i - 6) * v_t_unit(:, i)
            L_lin_others(:) = L_lin_others(:) + m_frag(i) * vtmp(:)
            call util_crossproduct(x_frag(:, i), vtmp(:), L(:))
            L_orb_others(:) = L_orb_others(:) + m_frag(i) * L(:)
         end if
      end do
      b(1:3) = -L_lin_others(:)
      b(4:6) = L_frag_orb(:) - L_orb_others(:)
      allocate(v_t_mag_output(nfrag))
      v_t_mag_output(1:6) = util_solve_linear_system(A, b, 6, lerr)
      if (present(v_t_mag_input)) v_t_mag_output(7:nfrag) = v_t_mag_input(:)

      return 
   end function solve_fragment_tangential_velocities

   subroutine set_fragment_radial_velocities(lerr)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy. This will minimize the difference between the fragment kinetic energy and the energy budget
      implicit none
      ! Arguments
      logical,                intent(out)   :: lerr
      ! Internals
      real(DP), parameter                   :: TOL = 1e-10_DP
      integer(I4B)                          :: i, j
      real(DP), dimension(:), allocatable   :: v_r_initial, v_r_sigma
      real(DP), dimension(:,:), allocatable :: v_r
      type(lambda_obj)                      :: objective_function

      ! Set the "target" ke_orbit_after (the value of the orbital kinetic energy that the fragments ought to have)
      
      allocate(v_r_initial, source=v_r_mag)
      ! Initialize radial velocity magnitudes with a random value that is approximately 10% of that found by distributing the kinetic energy equally
      allocate(v_r_sigma, source=v_r_mag)
      call random_number(v_r_sigma(1:nfrag))
      v_r_sigma(1:nfrag) = sqrt(1.0_DP + 2 * (v_r_sigma(1:nfrag) - 0.5_DP) * 1e-4_DP) 
      v_r_initial(1:nfrag) = v_r_sigma(1:nfrag) * sqrt(abs(ke_radial) / (2 * m_frag(1:nfrag) * nfrag)) 

      ! Initialize the lambda function using a structure constructor that calls the init method
      ! Minimize the ke objective function using the BFGS optimizer
      objective_function = lambda_obj(radial_objective_function)
      v_r_mag = util_minimize_bfgs(objective_function, nfrag, v_r_initial, TOL, lerr)
      ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
      vb_frag(:,1:nfrag) = vmag_to_vb(v_r_mag(1:nfrag), v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:)) 
      do i = 1, nfrag
         v_frag(:, i) = vb_frag(:, i) - vcom(:)
      end do
      ke_frag_orbit = 0.0_DP
      do i = 1, nfrag
         ke_frag_orbit = ke_frag_orbit + m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
      end do
      ke_frag_orbit = 0.5_DP * ke_frag_orbit
      !write(*,*) 'Radial'
      !write(*,*) 'Failure? ',lerr 
      !write(*,*) 'ke_frag_budget: ',ke_frag_budget
      !write(*,*) 'ke_frag_orbit : ',ke_frag_orbit
      !write(*,*) 'ke_frag_spin  : ',ke_frag_spin
      !write(*,*) 'ke_remainder  : ',ke_frag_budget - (ke_frag_orbit + ke_frag_spin)
      lerr = .false.

      return
   end subroutine set_fragment_radial_velocities

   function radial_objective_function(v_r_mag_input) result(fval) 
      !! Author: David A. Minton
      !!
      !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_r_mag_input   !! Unknown radial component of fragment velocity vector
      ! Result
      real(DP)                              :: fval      !! The objective function result, which is the square of the difference between the calculated fragment kinetic energy and our target
                                                         !! Minimizing this brings us closer to our objective
      ! Internals
      integer(I4B)                         :: i
      real(DP), dimension(:,:), allocatable :: v_shift

      allocate(v_shift, mold=vb_frag)
      v_shift(:,:) = vmag_to_vb(v_r_mag_input, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) 
      fval = 2 * ke_frag_budget 
      do i = 1, nfrag
         fval = fval - m_frag(i) * (Ip_frag(3, i) * rad_frag(i)**2 * dot_product(rot_frag(:, i), rot_frag(:, i)) + dot_product(v_shift(:, i), v_shift(:, i)))
      end do
      ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
      fval = (fval / (2 * ke_radial))**2

      return

   end function radial_objective_function

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
      implicit none
      integer(I4B) :: i
      real(DP), dimension(:), allocatable  :: m_frag_new, rad_frag_new
      real(DP), dimension(:,:), allocatable  :: xb_frag_new, vb_frag_new, Ip_frag_new, rot_frag_new
      real(DP) :: delta_r, delta_r_max
      real(DP), parameter :: ke_avg_deficit_target = 0.0_DP 

      ! Introduce a bit of noise in the radius determination so we don't just flip flop between similar failed positions
      call random_number(delta_r_max)
      delta_r_max = sum(radius(:)) * (1.0_DP + 2e-1_DP * (delta_r_max - 0.5_DP))
      if (try > 2) then
         ! Linearly interpolate the last two failed solution ke deficits to find a new distance value to try
         delta_r = (r_max_start - r_max_start_old) * (ke_avg_deficit_target - ke_avg_deficit_old) / (ke_avg_deficit - ke_avg_deficit_old)
         if (abs(delta_r) > delta_r_max) delta_r = sign(delta_r_max, delta_r)
      else
         delta_r = delta_r_max
      end if
      r_max_start_old = r_max_start
      r_max_start = r_max_start + delta_r ! The larger lever arm can help if the problem is in the angular momentum step
      if (f_spin > epsilon(1.0_DP)) then
         f_spin = f_spin / 2
      else
         f_spin = 0.0_DP
      end if
   end subroutine restructure_failed_fragments


end subroutine symba_frag_pos