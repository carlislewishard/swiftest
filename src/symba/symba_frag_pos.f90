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
   use module_interfaces, EXCEPT_THIS_ONE => symba_frag_pos
   implicit none
   ! Arguments
   type(user_input_parameters), intent(in)   :: param 
   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(:), intent(in)    :: family
   real(DP), intent(in)                      :: Qloss
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:), intent(in)        :: mass, radius, m_frag, rad_frag
   real(DP), dimension(:,:), intent(in)      :: Ip_frag
   real(DP), dimension(:,:), intent(out)     :: xb_frag, vb_frag, rot_frag
   logical, intent(out)                      :: lmerge ! Answers the question: Should this have been a merger instead?
   ! Internals
   real(DP), dimension(:,:), allocatable   :: x_frag, v_frag ! Fragment positions and velocities in the collision center of mass frame
   real(DP), dimension(NDIM, 2)            :: rot, L_orb 
   integer(I4B)                            :: i, j, nfrag, fam_size, istart
   real(DP), dimension(NDIM)               :: xcom, vcom, Ltot_before, Ltot_after, L_residual, L_spin_frag, L_target, L_frag
   real(DP)                                :: mtot, Lmag_before, Lmag_after
   real(DP)                                :: Etot_before, Etot_after, ke_before, pe_before
   real(DP)                                :: pe_after, ke_spin_before, ke_spin_after, ke_after, ke_family, ke_target, ke_frag
   real(DP), dimension(NDIM)               :: h, dx
   real(DP)                                :: rmag
   logical, dimension(:), allocatable      :: lexclude
   character(len=*), parameter             :: fmtlabel = "(A14,10(ES9.2,1X,:))"
   real(DP), dimension(NDIM)               :: x_cross_v, v_t_unit, h_unit, v_r_unit
   real(DP), dimension(:,:), allocatable   :: v_r, v_t
   
   associate(xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status, npl => symba_plA%helio%swiftest%nbody, name => symba_plA%helio%swiftest%id)

      allocate(x_frag, source=xb_frag)
      allocate(v_frag, source=vb_frag)
      allocate(v_r, mold=v_frag)
      allocate(v_t, mold=v_frag)
      fam_size = size(family)

      ! Find the center of mass of the collisional system	
      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 

      L_orb(:, :) = 0.0_DP
      ! Compute orbital angular momentum of pre-impact system
      do j = 1, 2
         call utiL_crossproduct(x(:, j) - xcom(:), v(:, j) - vcom(:), x_cross_v(:))
         L_orb(:, j) = mass(j) * x_cross_v(:)
      end do

      allocate(lexclude(npl))
      where (status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
         lexclude(1:npl) = .true.  
      elsewhere
         lexclude(1:npl) = .false. 
      end where

      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_before, ke_spin_before, pe_before, Ltot_before)
      Lmag_before = norm2(Ltot_before(:))
      Etot_before = ke_before + ke_spin_before + pe_before

      ! We need the original kinetic energy of just the pre-impact family members in order to balance the energy later
      ke_family = 0.0_DP
      do i = 1, fam_size
         ke_family = ke_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i))) !
         lexclude(family(i)) = .true. ! For all subsequent energy calculations the pre-impact family members will be replaced by the fragments
      end do
      ke_family = 0.5_DP * ke_family

      nfrag = size(m_frag)
      ! Initialize  positions and velocities of fragments that conserve angular momentum
      call symba_frag_pos_initialize_fragments(nfrag, xcom, vcom, x, v, L_orb, L_spin, mass, radius, m_frag, rad_frag, Ip_frag, x_frag, v_frag, rot_frag, lmerge)
      if (lmerge) return

      ! Energy calculation requires the fragments to be in the system barcyentric frame, so we need to temporarily shift them
      do i = 1, nfrag
         xb_frag(:,i) = x_frag(:,i) + xcom(:)
         vb_frag(:,i) = v_frag(:,i) + vcom(:)
      end do
     
      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot_after, &
         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
      Etot_after = ke_after + ke_spin_after + pe_after
      Lmag_after = norm2(Ltot_after(:))

      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,        "('              Energy normalized by |Etot_before|')")
      write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,        "('  First pass to get angular momentum ')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
                                         (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                         (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
                                         (pe_after - pe_before) / abs(Etot_before), &
                                         (Etot_after - Etot_before) / abs(Etot_before), &
                                         norm2(Ltot_after - Ltot_before) / Lmag_before
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,        "('  Second pass to get energy ')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' Q_loss      |',-Qloss / abs(Etot_before)
      write(*,        "(' ---------------------------------------------------------------------------')")

      
      ! Set the "target" ke_after (the value of the orbital kinetic energy that the fragments ought to have)
      ke_target = ke_family + (ke_spin_before - ke_spin_after) + (pe_before - pe_after) - Qloss
      L_target(:) = 0.0_DP
      do i = 1, nfrag
         call util_crossproduct(x_frag(:, i), v_frag(:, i), L_frag(:))
         L_target(:) = L_target(:) + L_frag(:)
      end do
      call symba_frag_pos_kinetic_energy(m_frag, ke_target, L_target, x_frag, v_frag, lmerge)
      
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' T_family    |',ke_family / abs(Etot_before)
      write(*,fmtlabel) ' T_frag targ |',ke_target / abs(Etot_before)

      ! Shift the fragments into the system barycenter frame
      do i = 1, nfrag
         xb_frag(:,i) = x_frag(:, i) + xcom(:)
         vb_frag(:,i) = v_frag(:, i) + vcom(:)
      end do

      ke_frag = 0._DP
      do i = 1, nfrag
         ke_frag = ke_frag + 0.5_DP * m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
      end do
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' T_frag calc |',ke_frag / abs(Etot_before)
      write(*,fmtlabel) ' residual    |',1.0_DP - ke_frag / ke_target

      ! Calculate the new energy of the system of fragments
      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot_after,&
            nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
      Etot_after = ke_after + ke_spin_after + pe_after
      Lmag_after = norm2(Ltot_after(:))

      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
                                        (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                        (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
                                        (pe_after - pe_before) / abs(Etot_before), &
                                        (Etot_after - Etot_before) / abs(Etot_before), &
                                        norm2(Ltot_after - Ltot_before) / Lmag_before
   
      lmerge = lmerge .or. ((Etot_after - Etot_before) / abs(Etot_before) > 0._DP) 
      !if (.not.lmerge) then
      !   L_residual(:) = Ltot_before(:) - Ltot_after(:)
      !   L_spin_frag(:) = L_residual(:) / nfrag
      !   do i = 1, nfrag
      !      rot_frag(:,i) = rot_frag(:,i) + L_spin_frag(:) / (Ip_frag(:, i) * m_frag(i) * rad_frag(i)**2)
      !   end do
      !end if

!      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot_after,&
!         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
!         Etot_after = ke_after + ke_spin_after + pe_after
!      L_residual(:) = Ltot_before(:) - Ltot_after(:)
!      Lmag_after = norm2(Ltot_after(:))
!
!      write(*,        "(' ---------------------------------------------------------------------------')")
!      write(*,        "('  Third pass for correcting any residual angular momentum ')")
!      write(*,        "(' ---------------------------------------------------------------------------')")
!      !write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
!      !write(*,        "(' ---------------------------------------------------------------------------')")
!      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
!                                       (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
!                                       (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
!                                       (pe_after - pe_before) / abs(Etot_before), &
!                                       (Etot_after - Etot_before) / abs(Etot_before), &
!                                       norm2(Ltot_after - Ltot_before) / Lmag_before
!      write(*,        "(' ---------------------------------------------------------------------------')")
      !****************************************************************************************************************

   end associate
   return 

   contains

   subroutine symba_frag_pos_initialize_fragments(nfrag, xcom, vcom, x, v, L_orb, L_spin, mass, radius, m_frag, rad_frag, Ip_frag, x_frag, v_frag, rot_frag, lmerge)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
      !! pre-impact angular momentum. They are distributed on an ellipse surrounding the center of mass.
      !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.
      implicit none
      ! Arguments
      integer(I4B),             intent(in)    :: nfrag                    !! Number of collisional fragments
      real(DP), dimension(:),   intent(in)    :: xcom, vcom               !! Center of mass position and velocity in the system barycenter frame
      real(DP), dimension(:,:), intent(in)    :: x, v, L_orb, L_spin     !! Pre-impact angular momentum vectors
      real(DP), dimension(:),   intent(in)    :: mass, radius             !! Pre-impact masses and radii
      real(DP), dimension(:),   intent(in)    :: m_frag, rad_frag         !! Fragment masses and radii
      real(DP), dimension(:,:), intent(in)    :: Ip_frag                  !! Fragment prinicpal moments of inertia
      real(DP), dimension(:,:), intent(out)   :: x_frag, v_frag, rot_frag !! Fragment position, velocities, and spin states
      logical, intent(out)                    :: lmerge
      ! Internals
      real(DP)                                :: mtot, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
      real(DP), dimension(NDIM)               :: Ltot, delta_r, delta_v
      real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
      integer(I4B)                            :: i

      ! Now create the fragment distribution
      mtot = sum(mass(:))

      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      Ltot = L_spin(:,1) + L_spin(:,2) + L_orb(:, 1) + L_orb(:, 2)

      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))     
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))

      ! We will initialize fragments on a plane defined by the pre-impact system, with the y-axis aligned with the angular momentum vector
      ! and the z-axis aligned with the pre-impact distance vector.
      y_col_unit = Ltot(:) / norm2(Ltot(:))
      z_col_unit(:) = delta_r(:) / r_col_norm  
      ! The cross product of the z- by x-axis will give us the y-axis
      call util_crossproduct(y_col_unit, z_col_unit, x_col_unit)

      ! The angular spacing of fragments on the ellipse - We will only use half the ellipse
      theta = 2 * PI / nfrag

      ! Re-normalize position and velocity vectors by the fragment number so that for our initial guess we weight each
      ! fragment position by the mass and assume equipartition of energy for the velocity
      r_col_norm = 1.5_DP * 2 * sum(rad_frag(:)) / theta / nfrag 

      ! We will treat the first fragment of the list as a special case.
      x_frag(:, 1) = -z_col_unit(:) 
      call random_number(x_frag(:,2:nfrag)) 
      
      x_frag(:, :) = x_frag(:, :) * r_col_norm 
      call symba_frag_pos_com_adjust(xcom, m_frag, x_frag)
      v_frag(:,:) = 0._DP

      call symba_frag_pos_ang_mtm(nfrag, xcom, vcom, L_orb, L_spin, m_frag, rad_frag, Ip_frag, x_frag, v_frag, rot_frag, lmerge)

      return
   end subroutine symba_frag_pos_initialize_fragments

   subroutine symba_frag_pos_ang_mtm(nfrag, xcom, vcom, L_orb, L_spin, m_frag, rad_frag, Ip_frag, x_frag, v_frag, rot_frag, lmerge)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
      implicit none
      ! Arguments
      integer(I4B),             intent(in)    :: nfrag            !! Number of collisional fragments
      real(DP), dimension(:),   intent(in)    :: xcom, vcom       !! Center of mass position and velocity in the system barycenter frame
      real(DP), dimension(:,:), intent(in)    :: L_orb, L_spin     !! Pre-impact position, velocity, and spin states, Ip_frag
      real(DP), dimension(:),   intent(in)    :: m_frag, rad_frag !! Fragment masses and radii
      real(DP), dimension(:,:), intent(in)    :: Ip_frag, x_frag  !! Fragment prinicpal moments of inertia and position vectors
      real(DP), dimension(:,:), intent(out)   :: v_frag, rot_frag !! Fragment velocities and spin states
      logical, intent(out)                    :: lmerge
      ! Internals
      integer(I4B)                            :: i
      real(DP)                                :: L_orb_mag
      real(DP), dimension(NDIM)               :: L_orb_old, L_spin_frag, L_spin_new
      real(DP), parameter                     :: f_spin = 0.20_DP !! Fraction of pre-impact orbital angular momentum that is converted to fragment spin
      real(DP), dimension(2 * NDIM, 2 * NDIM) :: A ! LHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(2 * NDIM)           :: b  ! RHS of linear equation used to solve for momentum constraint in Gauss elimination code
      real(DP), dimension(:,:), allocatable   :: v_t_unit, r_unit
      real(DP), dimension(:), allocatable     :: v_t_mag, rmag
      real(DP), dimension(NDIM)               :: L_lin_others, L_orb_others, L_unit, L_frag, L_sigma

      allocate(v_t_mag, mold=m_frag)
      allocate(rmag, mold=m_frag)
      allocate(r_unit, mold=v_frag)
      allocate(v_t_unit, mold=v_frag)
      
      L_orb_old(:) = L_orb(:, 1) + L_orb(:, 2)
   
      ! Divide up the pre-impact spin angular momentum equally between the various bodies by mass
      L_spin_new(:) = L_spin(:,1) + L_spin(:, 2) + f_spin * L_orb_old(:)
      L_spin_frag(:) = L_spin_new(:) / nfrag
      do i = 1, nfrag
         rot_frag(:,i) = L_spin_frag(:) / (Ip_frag(:, i) * m_frag(i) * rad_frag(i)**2)
      end do
      L_orb_old(:) = L_orb_old(:) * (1.0_DP - f_spin)
      L_unit(:) = L_orb_old(:) / norm2(L_orb_old(:))

      L_orb_mag = norm2(L_orb_old(:)) 
      ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
      ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
      ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
      L_lin_others(:) = 0.0_DP
      L_orb_others(:) = 0.0_DP
      do i = 1, nfrag
         rmag(i) = norm2(x_frag(:, i))
         r_unit(:, i) = x_frag(:, i) / rmag(i)
         call random_number(L_sigma(:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                        ! otherwise we can get an ill-conditioned system
         L_frag(:) = L_unit(:) + 2e-3_DP * (L_sigma(:) - 0.5_DP)
         L_frag(:) = L_frag(:) / norm2(L_frag(:)) 
         call util_crossproduct(L_frag(:), r_unit(:, i), v_t_unit(:, i))
         if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
            A(1:3, i) = m_frag(i) * v_t_unit(:, i) 
            call util_crossproduct(r_unit(:, i), v_t_unit(:, i), L_frag(:))
            A(4:6, i) = m_frag(i) * rmag(i) * L_frag(:)
         else !For the remining bodies, distribute the angular momentum equally amongs them
            v_t_mag(i) = L_orb_mag / (m_frag(i) * rmag(i) * nfrag)
            v_frag(:, i) = v_t_mag(i) * v_t_unit(:, i)
            L_lin_others(:) = L_lin_others(:) + m_frag(i) * v_frag(:, i)
            call util_crossproduct(x_frag(:, i), v_frag(:, i), L_frag(:))
            L_orb_others(:) = L_orb_others(:) + m_frag(i) * L_frag(:)
         end if
      end do
      b(1:3) = -L_lin_others(:)
      b(4:6) = L_orb_old(:) - L_orb_others(:)
      v_t_mag(1:6) = util_solve_linear_system(A, b, 6, lmerge)
      if (lmerge) return
      do i = 1, 6
         v_frag(:, i) = v_t_mag(i) * v_t_unit(:, i)
      end do

      return 
   end subroutine symba_frag_pos_ang_mtm

   subroutine symba_frag_pos_kinetic_energy(m_frag, ke_target, L_target, x_frag, v_frag, lmerge)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy.
      !! It will check that we don't end up with negative energy (bound system). If so, we'll set the fragment velocities to
      !! zero in the center of mass frame and indicate the the fragmentation should instead by a merger.
      !! It takes in the initial "guess" of velocities and solve for the a scaling factor applied to the radial component wrt the
      !! center of mass frame needed to correct the kinetic energy of the fragments in the system barycenter frame to match that of 
      !! the target kinetic energy required to satisfy the constraints.
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: m_frag     !! Fragment masses
      real(DP), dimension(:,:), intent(in)    :: x_frag     !! Fragment position vectors
      real(DP),                 intent(in)    :: ke_target  !! Target kinetic energy 
      real(DP), dimension(:),   intent(in)    :: L_target   !! Target kinetic energy 
      real(DP), dimension(:,:), intent(inout) :: v_frag      !! Fragment position and velocities in the center of mass frame   
      logical,                  intent(out)   :: lmerge

      ! Internals
      real(DP)                              :: mtot           !! Total mass of fragments
      integer(I4B)                          :: i, nfrag, neval
      real(DP), parameter                   :: TOL = 1e-9_DP
      real(DP), dimension(:), allocatable   :: vflat
      logical                               :: lerr

         
      nfrag = size(m_frag)
      mtot = sum(m_frag)

      ! Initialize the lambda function using a structure constructor that calls the init method
      ! Minimize error using the BFGS optimizer
      vflat = util_minimize_bfgs(symba_frag_lambda(lambda=symba_frag_pos_ke_objective_function, v_frag=v_frag, x_frag=x_frag, m_frag=m_frag, L_target=L_target, ke_target=ke_target), &
                                 NDIM*nfrag, reshape(v_frag,[NDIM*nfrag]), TOL, lerr)
      if (lerr) then
         ! No solution exists for this case, so we need to indicate that this should be a merge
         ! This may happen due to setting the tangential velocities too high when setting the angular momentum constraint
         lmerge = .true.
         v_frag(:,:) = 0.0_DP
      else
         v_frag(:,:) = reshape(vflat,shape(v_frag))
      end if

      return
   end subroutine symba_frag_pos_kinetic_energy

   function symba_frag_pos_ke_objective_function(vflat, v_frag, x_frag, m_frag, L_target, ke_target) result(fnorm)
      ! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in) :: vflat          !! Unrolled unknown velocity vectors
      real(DP), dimension(:,:), intent(in) :: v_frag, x_frag !! Velocity and position vectors
      real(DP), dimension(:),   intent(in) :: m_frag         !! Fragment masses
      real(DP), dimension(:),   intent(in) :: L_target       !! Target orbital momentum
      real(DP),                 intent(in) :: ke_target      !! Target kinetic energ
      ! Result
      real(DP)                             :: fnorm     !! The objective function result: norm of the vector composed of the tangential momentum and energy
                                                        !! Minimizing this brings us closer to our objective
      ! Internals
      real(DP), dimension(7)  :: f_vec
      integer(I4B) :: i, nfrag, nsol
      real(DP), dimension(NDIM) :: L
      real(DP), dimension(:,:), allocatable :: v_sol

      nfrag = size(m_frag)
      nsol = size(vflat) / NDIM
      allocate(v_sol(NDIM,nsol))
      v_sol = reshape(vflat, shape(v_sol))

      ! Linear momentum constraint
      f_vec(1) = sum(m_frag(:) * [v_sol(1, 1:nsol), v_frag(1, nsol + 1:nfrag)])
      f_vec(2) = sum(m_frag(:) * [v_sol(2, 1:nsol), v_frag(2, nsol + 1:nfrag)])
      f_vec(3) = sum(m_frag(:) * [v_sol(3, 1:nsol), v_frag(3, nsol + 1:nfrag)])
      ! Angular momentum and kinetic energy constraints
      f_vec(4:6) = -L_target(:)
      f_vec(7) = -ke_target
      do i = 1, nsol
         call util_crossproduct(x_frag(:,i), v_sol(:, i), L(:))
         f_vec(4:6) = f_vec(4:6) + m_frag(i) * L(:)
         f_vec(7)   = f_vec(7)   + 0.5_DP * m_frag(i) * dot_product(v_sol(:,i), v_sol(:, i))
      end do
      do i = nsol + 1, nfrag
         call util_crossproduct(x_frag(:,i), v_frag(:, i), L(:))
         f_vec(4:6) = f_vec(4:6) + m_frag(i) * L(:)
         f_vec(7)   = f_vec(7) + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:, i))
      end do

      fnorm = norm2(f_vec(:)) 

      return

   end function symba_frag_pos_ke_objective_function

   subroutine symba_frag_pos_com_adjust(vec_com, m_frag, vec_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the position or velocity of the fragments as needed to align them with the original trajectory center of mass.
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: vec_com   !! Center of mass position or velocity in the system barycenter frame
      real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame

      ! Internals
      real(DP), dimension(NDIM)               :: mvec_frag, COM_offset
      real(DP)                                :: mtot
      integer(I4B)                            :: i, nfrag
         
      mtot = sum(m_frag(:))
      nfrag = size(m_frag(:))
      mvec_frag(:) = 0.0_DP

      do i = 1, nfrag
         mvec_frag = mvec_frag(:) + (vec_frag(:,i) + vec_com(:)) * m_frag(i)
      end do
      COM_offset(:) = vec_com(:) - mvec_frag(:) / mtot
      do i = 1, nfrag 
         vec_frag(:, i) = vec_frag(:, i) + COM_offset(:)
      end do

      return
   end subroutine symba_frag_pos_com_adjust

   subroutine symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_orbit, ke_spin, pe, Ltot, nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy, including all bodies in the symba_plA list that do not have a corresponding value of the lexclude array that is true
      !! and optionally including fragments.
      use module_swiftestalloc
      implicit none
      ! Arguments
      integer(I4B), intent(inout) :: npl
      type(symba_pl), intent(inout) :: symba_plA
      logical, dimension(:), allocatable, intent(inout) :: lexclude
      real(DP), intent(out) :: ke_orbit, ke_spin, pe
      real(DP), dimension(:), intent(out)      :: Ltot
      integer(I4B), intent(in), optional :: nfrag
      real(DP), dimension(:), intent(in), optional :: m_frag, rad_frag
      real(DP), dimension(:,:), intent(in), optional :: Ip_frag, xb_frag, vb_frag, rot_frag
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
      lk_plpl = allocated(symba_plA%helio%swiftest%k_plpl)
      if (lk_plpl) deallocate(symba_plA%helio%swiftest%k_plpl) 
      if (present(nfrag)) then ! Temporarily expand the planet list to feed it into symba_energy
         npl_new = npl + nfrag
      else
         npl_new  = npl
      end if
      call symba_pl_allocate(symba_plwksp, npl_new)

      ! Copy over old data
      symba_plwksp%helio%swiftest%id(1:npl) = symba_plA%helio%swiftest%id(1:npl)
      symba_plwksp%helio%swiftest%status(1:npl) = symba_plA%helio%swiftest%status(1:npl)
      symba_plwksp%helio%swiftest%mass(1:npl) = symba_plA%helio%swiftest%mass(1:npl)
      symba_plwksp%helio%swiftest%radius(1:npl) = symba_plA%helio%swiftest%radius(1:npl)
      symba_plwksp%helio%swiftest%xh(:,1:npl) = symba_plA%helio%swiftest%xh(:,1:npl)
      symba_plwksp%helio%swiftest%vh(:,1:npl) = symba_plA%helio%swiftest%vh(:,1:npl)
      symba_plwksp%helio%swiftest%rhill(1:npl) = symba_plA%helio%swiftest%rhill(1:npl)
      symba_plwksp%helio%swiftest%xb(:,1:npl) = symba_plA%helio%swiftest%xb(:,1:npl)
      symba_plwksp%helio%swiftest%vb(:,1:npl) = symba_plA%helio%swiftest%vb(:,1:npl)
      symba_plwksp%helio%swiftest%rot(:,1:npl) = symba_plA%helio%swiftest%rot(:,1:npl)
      symba_plwksp%helio%swiftest%Ip(:,1:npl) = symba_plA%helio%swiftest%Ip(:,1:npl)

      if (present(nfrag)) then ! Append the fragments if they are included
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
      end if

      where (lexclude(1:npl))
         symba_plwksp%helio%swiftest%status(1:npl) = INACTIVE
      end where

      nplm = count(symba_plwksp%helio%swiftest%mass > param%mtiny)
      call util_dist_index_plpl(npl_new, nplm, symba_plwksp)
      call symba_energy_eucl(npl_new, symba_plwksp, param%j2rp2, param%j4rp4, ke_orbit, ke_spin, pe, te, Ltot)

      ! Restore the big array
      deallocate(symba_plwksp%helio%swiftest%k_plpl) 
      nplm = count(symba_plA%helio%swiftest%mass > param%mtiny)
      if (lk_plpl) call util_dist_index_plpl(npl, nplm, symba_plA)

      return
   end subroutine symba_frag_pos_energy_calc
 




end subroutine symba_frag_pos
