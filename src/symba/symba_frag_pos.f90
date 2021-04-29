subroutine symba_frag_pos (param, symba_plA, idx_parents, x, v, L_spin, Ip, mass, radius, &
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
   integer(I4B), dimension(:), intent(in)    :: idx_parents
   real(DP), intent(in)                      :: Qloss
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:), intent(in)        :: mass, radius, m_frag, rad_frag
   real(DP), dimension(:,:), intent(in)      :: Ip_frag
   real(DP), dimension(:,:), intent(out)     :: xb_frag, vb_frag, rot_frag
   logical, intent(out)                      :: lmerge ! Answers the question: Should this have been a merger instead?
   ! Internals
   real(DP), dimension(:,:), allocatable   :: x_frag, v_frag ! Fragment positions and velocities in the collision center of mass frame
   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, nfrag, fam_size, istart, npl
   real(DP), dimension(NDIM)               :: xcom, vcom, Ltot
   real(DP)                                :: mtot, Ltot_before, Ltot_after
   real(DP)                                :: Etot_before, Etot_after, ke_before, pe_before
   real(DP)                                :: pe_after, ke_spin_before, ke_spin_after, ke_after, ke_family, ke_target
   real(DP), dimension(NDIM)               :: h, dx
   integer(I4B), dimension(:), allocatable :: family
   real(DP)                                :: rmag
   logical, dimension(:), allocatable      :: lfamily
   logical, dimension(:), allocatable      :: lexclude
   character(len=*), parameter             :: fmtlabel = "(A14,6(ES9.2,1X,:))"
   
   associate(nchild1 => symba_plA%kin(idx_parents(1))%nchild, nchild2 => symba_plA%kin(idx_parents(2))%nchild, &
             xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status, npl => symba_plA%helio%swiftest%nbody)

      allocate(x_frag, source=xb_frag)
      allocate(v_frag, source=vb_frag)
      allocate(lexclude(npl))
      lexclude(:) = status(1:npl) == INACTIVE
      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_before, ke_spin_before, pe_before, Ltot)
      Ltot_before = norm2(Ltot(:))
      Etot_before = ke_before + ke_spin_before + pe_before

      ! Find the center of mass of the collisional system
      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 

      ! Make the list of family members (bodies involved in the collision)
      fam_size = 2 + nchild1 + nchild2
      allocate(family(fam_size))
      allocate(lfamily(fam_size))
      family(1) = idx_parents(1)
      family(2) = idx_parents(2)
      istart = 2 + nchild1

      if (nchild1 > 0) family(3:istart) = symba_plA%kin(idx_parents(1))%child(1:nchild1)
      if (nchild2 > 0) family(istart+1:istart+1+nchild2) = symba_plA%kin(idx_parents(2))%child(1:nchild2)
      lfamily(:) = (status(family(:)) == ACTIVE) .or. (status(family(:)) == COLLISION) 
      fam_size = count(lfamily(:))
      family(:) = pack(family(:), lfamily(:))
      ! We need the original kinetic energy of just the pre-impact family members in order to balance the energy later
      ke_family = 0.0_DP
      do i = 1, fam_size
         ke_family = ke_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i))) ! 
         lexclude(family(i)) = .true. ! Exclude pre-impact family members from subsequent energy calculations
      end do
      ke_family = 0.5_DP * ke_family

      call symba_frag_pos_initialize_fragments(xcom, vcom, x, v, L_spin, mass, m_frag, x_frag, v_frag)   

      ! Conserve the system angular momentum (but not energy)
      call symba_frag_pos_ang_mtm(xcom, vcom, x, v, mass, radius, L_spin, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)

      ! Energy calculation requires the fragments to be in the system barcyentric frame, so we need to temporarily shift them
      nfrag = size(m_frag)
      do i = 1, nfrag
         xb_frag(:,i) = x_frag(:,i) + xcom(:)
         vb_frag(:,i) = v_frag(:,i) + vcom(:)
      end do
     
      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot, &
         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
      Etot_after = ke_after + ke_spin_after + pe_after
      Ltot_after = norm2(Ltot(:))

      write(*,        "('              Energy normalized by |Etot_before|')")
      write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' original    |',ke_before / abs(Etot_before), &
                                         ke_spin_before / abs(Etot_before), &
                                         (ke_before + ke_spin_before) / abs(Etot_before), &
                                         pe_before / abs(Etot_before), &
                                         Etot_before / abs(Etot_before), &
                                         Ltot_before / Ltot_before
      write(*,fmtlabel) ' first pass  |',ke_after / abs(Etot_before), &
                                         ke_spin_after / abs(Etot_before), &
                                         (ke_after + ke_spin_after) / abs(Etot_before), &
                                         pe_after / abs(Etot_before), &
                                         Etot_after / abs(Etot_before), &
                                         Ltot_after / Ltot_before
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
                                         (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                         (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
                                         (pe_after - pe_before) / abs(Etot_before), &
                                         (Etot_after - Etot_before) / abs(Etot_before), &
                                         (Ltot_after - Ltot_before) / Ltot_before
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' Q_loss      |',-Qloss / abs(Etot_before)
      write(*,        "(' ---------------------------------------------------------------------------')")

      ! Set the "target" ke_after (the value of the orbital kinetic energy that the fragments ought to have)
      ke_target = ke_family + (ke_spin_before - ke_spin_after) + (pe_before - pe_after) - Qloss
      call symba_frag_pos_kinetic_energy(xcom, vcom, m_frag, x_frag, v_frag, ke_target, lmerge)
      write(*,fmtlabel) ' target      |',ke_target / abs(Etot_before)
      write(*,        "(' ---------------------------------------------------------------------------')")

      ! Shift the fragments into the system barycenter frame
      do i = 1, nfrag
         xb_frag(:,i) = x_frag(:, i) + xcom(:)
         vb_frag(:,i) = v_frag(:, i) + vcom(:)
      end do

      ! REMOVE THE FOLLOWING AFTER TESTING
      !****************************************************************************************************************
      ! Calculate the new energy of the system of fragments
      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot, &
         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
      Etot_after = ke_after + ke_spin_after + pe_after
      Ltot_after = norm2(Ltot(:))
     
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' final       |',ke_after / abs(Etot_before), &
                                         ke_spin_after / abs(Etot_before), &
                                         (ke_after + ke_spin_after) / abs(Etot_before), &
                                         pe_after / abs(Etot_before), &
                                         Etot_after / abs(Etot_before), &
                                         Ltot_after / Ltot_before
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,        "(' (final T_orb should be the same as target T_orb)')")
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
                                         (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
                                         (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
                                         (pe_after - pe_before) / abs(Etot_before), &
                                         (Etot_after - Etot_before) / abs(Etot_before), &
                                         (Ltot_after - Ltot_before) / Ltot_before
      write(*,        "(' ---------------------------------------------------------------------------')")
      write(*,*)   
      !****************************************************************************************************************

   end associate
   return 

   contains

   subroutine symba_frag_pos_initialize_fragments(xcom, vcom, x, v, L_spin, mass, m_frag, x_frag, v_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
      !! pre-impact angular momentum. They are distributed on a circle surrounding the center of mass and with velocities pointing outward.
      !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: xcom, vcom      !! Center of mass position and velocity in the system barycenter frame
      real(DP), dimension(:,:), intent(in)    :: x, v, L_spin    !! Pre-impact spins
      real(DP), dimension(:), intent(in)      :: mass            !! Pre-impact masses
      real(DP), dimension(:), intent(in)      :: m_frag          !! Fragment masses
      real(DP), dimension(:,:), intent(out)   :: x_frag, v_frag  !! Fragment position and velocities

      ! Internals
      integer(I4B), save                      :: thetashift = 0
      integer(I4B), parameter                 :: SHIFTMAX = 9
      real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
      real(DP), dimension(NDIM)               :: Ltot, xc, vc, x_cross_v, delta_r, delta_v
      real(DP), dimension(NDIM)               :: r_col_unit_vec, v_col_unit_vec, v_plane_unit_vec
      integer(I4B)                            :: i, nfrag

      ! Now create the fragment distribution
      nfrag = size(m_frag(:))
      mtot = sum(mass(:))

      ! Calculate the position of each fragment 
      ! Theta is a phase shift value that ensures that successive nearby collisions in a single step are rotated to avoid possible overlap
      theta = (2 * PI) / nfrag
      ! Shifts the starting circle of fragments around so that multiple fragments generated 
      ! from a single collision in a single time step don't pile up on top of each other
      phase_ang = theta * thetashift / SHIFTMAX
      thetashift = thetashift + 1
      IF (thetashift >= shiftmax) thetashift = 0

      ! Theta is a phase shift value that ensures that successive nearby collisions in a single step are rotated to avoid possible overlap
      theta = (2 * PI) / nfrag
      ! Shifts the starting circle of fragments around so that multiple fragments generated 
      ! from a single collision in a single time step don't pile up on top of each other
      phase_ang = theta * thetashift / SHIFTMAX
      thetashift = thetashift + 1
      IF (thetashift >= shiftmax) thetashift = 0

      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      Ltot = L_spin(:,1) + L_spin(:,2)
      do i = 1, 2
         xc(:) = x(:, i) - xcom(:)
         vc(:) = v(:, i) - vcom(:)
         call utiL_crossproduct(xc(:), vc(:), x_cross_v(:))
         Ltot(:) = Ltot(:) + mass(i) * x_cross_v(:)
      end do

      ! Relative position and velocity vectors of the two impacting "clouds" 
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))
      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))               

      ! Calculate the triple product to get the plane of the fragment distribution
      call util_crossproduct(Ltot,delta_v,v_plane_unit_vec)
      v_plane_unit_vec(:) = v_plane_unit_vec(:) / norm2(v_plane_unit_vec(:))

      v_col_unit_vec(:) = delta_v(:) / v_col_norm 
      r_col_unit_vec(:) = delta_r(:) / norm2(delta_r(:)) ! unit vector of collision distance

      ! Re-normalize position and velocity vectors by the fragment number so that for our initial guess we weight each
      ! fragment position by the mass and assume equipartition of energy for the velocity
      r_col_norm = max(2 * r_col_norm, 2 * sum(radius(:))) / nfrag ! To ensure that the new fragments aren't overlapping we will pick an initial starting radius 
                                                                   ! that is the bigger of: 2x the initial separation or 2x the mutual radius. 
      v_col_norm = v_col_norm / sqrt(1.0_DP * nfrag)
      do i = 1, nfrag
         ! Place the fragments on the collision plane at a distance proportional to mass wrt the collisional barycenter
         ! This gets updated later after the new potential energy is calculated
         r_frag_norm = r_col_norm * mtot / m_frag(i) 

         x_frag(:,i) =  r_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)  + &
                                       (sin(phase_ang + theta * i)) * v_plane_unit_vec(:)) 
                        
         ! Apply a simple mass weighting first to ensure that the velocity follows the barycenter
         ! This gets updated later after the new potential and kinetic energy is calcualted
         v_frag_norm = v_col_norm * sqrt(mtot / m_frag(i))
         v_frag(:,i) =  v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:) + &
                                       (sin(phase_ang + theta * i)) * v_plane_unit_vec(:)) 
      end do

      return
   end subroutine symba_frag_pos_initialize_fragments

   subroutine symba_frag_pos_ang_mtm(xcom, vcom, x, v, mass, radius, L_spin, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
      implicit none
   
      real(DP), dimension(:),   intent(in)    :: xcom, vcom, mass, radius
      real(DP), dimension(:,:), intent(in)    :: x, v, L_spin
      real(DP), dimension(:,:), intent(in)    :: Ip_frag
      real(DP), dimension(:),   intent(in)    :: m_frag, rad_frag
      real(DP), dimension(:,:), intent(inout) :: x_frag, v_frag, rot_frag
   
      real(DP), dimension(NDIM, 2)            :: rot, Ip
      integer(I4B)                            :: i, j, nfrag
      real(DP)                                :: mtot
      real(DP), dimension(NDIM)               :: xc, vc, x_cross_v, v_phi_unit, h_unit, v_r_unit
      real(DP), dimension(NDIM)               :: L_orb_old, L_spin_old, L_orb_new, L_spin_frag, L_residual, L_spin_new
   
      nfrag = size(m_frag)
      mtot = sum(mass(:))
      
      L_spin_old(:) = L_spin(:,1) + L_spin(:,2)
      L_orb_old(:) = 0.0_DP
      ! Compute orbital angular momentum of pre-impact system
      do j = 1, 2
         xc(:) = x(:, j) - xcom(:)
         vc(:) = v(:, j) - vcom(:)
         call utiL_crossproduct(xc(:), vc(:), x_cross_v(:))
         L_orb_old(:) = L_orb_old(:) + mass(j) * x_cross_v(:)
      end do
   
      ! Divide up the pre-impact spin angular momentum equally between the various bodies by mass
      L_spin_new(:) = L_spin_old(:)
      do i = 1, nfrag
         L_spin_frag(:) = L_spin_new(:) * m_frag(i) / mtot
         rot_frag(:,i) = L_spin_frag(:) / (Ip_frag(3, i) * m_frag(i) * rad_frag(i)**2)
      end do
      
      do i = 1, nfrag
         L_residual(:) = L_orb_old(:) * m_frag(i) / mtot
         h_unit(:) = L_orb_old(:) / norm2(L_orb_old(:))
         v_r_unit(:) = x_frag(:,i) / norm2(x_frag(:,i))
         call util_crossproduct(h_unit(:), v_r_unit(:), v_phi_unit(:))  ! make a unit vector in the tangential velocity direction
         v_frag(:,i) = v_frag(:,i) + norm2(L_residual(:)) / m_frag(i) / norm2(x_frag(:,i)) * v_phi_unit(:) ! Distribute the angular momentum equally amongst the fragments
      end do

      call symba_frag_pos_com_adjust(xcom, vcom, m_frag, x_frag, v_frag)
   
      return 
   end subroutine symba_frag_pos_ang_mtm

   subroutine symba_frag_pos_kinetic_energy(xcom, vcom, m_frag, x_frag, v_frag, ke_target, lmerge)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinent 
      !! It will check that we don't end up with negative energy (bound system). If so, we'll set the fragment velocities to
      !! zero in the center of mass frame and indicate the the fragmentation should instead by a merger.
      !! It takes in the initial "guess" of velocities and solve for the a scaling factor applied to the radial component wrt the
      !! center of mass frame needed to correct the kinetic energy of the fragments in the system barycenter frame to match that of 
      !! the target kinetic energy required to satisfy the constraints.
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: xcom, vcom      !! Center of mass position and velocity in the system barycenter frame
      real(DP), dimension(:),   intent(in)    :: m_frag          !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: x_frag, v_frag  !! Fragment position and velocities in the center of mass frame   
      real(DP), intent(in)                    :: ke_target        !! Target kinetic energy 
      logical, intent(out)                    :: lmerge

      ! Internals
      real(DP)                                :: f_corrected, A, B, C, rterm, mtot
      integer(I4B)                            :: i, nfrag
      real(DP), dimension(:,:), allocatable   :: v_r, v_phi
      real(DP), dimension(NDIM)               :: x_cross_v, v_phi_unit, h_unit, v_r_unit
         
      nfrag = size(m_frag)

      allocate(v_r(NDIM,nfrag))
      allocate(v_phi(NDIM,nfrag))

      do i = 1, nfrag
         call utiL_crossproduct(x_frag(:,i), v_frag(:,i), x_cross_v(:))
         h_unit(:) = x_cross_v(:) / norm2(x_cross_v(:))
         v_r_unit(:) = x_frag(:,i) / norm2(x_frag(:, i))
         call utiL_crossproduct(h_unit(:), v_r_unit(:), v_phi_unit(:))
         v_r(:,i) = dot_product(v_frag(:,i), v_r_unit(:)) * v_r_unit(:)
         v_phi(:,i) = dot_product(v_frag(:,i), v_phi_unit(:)) * v_phi_unit(:)
      end do

      A = 0.0_DP
      B = 0.0_DP
      C = 0.0_DP
      do i = 1, nfrag
         A = A + m_frag(i) * dot_product(v_r(:,i), v_r(:,i))
         B = B + m_frag(i) * dot_product(v_r(:,i), vcom(:))
         C = C + m_frag(i) * (0.5_DP * (dot_product(v_phi(:,i), v_phi(:,i)) + dot_product(vcom(:), vcom(:))) + dot_product(v_phi(:,i), vcom(:)))
      end do
      A = 0.5_DP * A
      C = C - ke_target
      rterm = B**2 - 4 * A * C
      if (rterm > 0.0_DP) then
         f_corrected = (-B + sqrt(rterm)) / (2 * A)
         lmerge = .false.
      else
         f_corrected = 0.0_DP
         lmerge = .true.
      end if

      ! Shift the fragments into the system barycenter frame
      v_frag(:,:) = f_corrected * v_r(:, :) + v_phi(:, :) 

      call symba_frag_pos_com_adjust(xcom, vcom, m_frag, x_frag, v_frag)

      write(*,fmtlabel) ' f_corrected |',f_corrected

      return
   end subroutine symba_frag_pos_kinetic_energy

   subroutine symba_frag_pos_com_adjust(xcom, vcom, m_frag, x_frag, v_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjust the position and velocity of the fragments as needed to align them with the original trajectory center of mass
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: xcom, vcom      !! Center of mass position and velocity in the system barycenter frame
      real(DP), dimension(:),   intent(in)    :: m_frag          !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: x_frag, v_frag  !! Fragment position and velocities in the center of mass frame

      ! Internals
      real(DP), dimension(NDIM)               :: mx_frag, mv_frag, COM_offset_x, COM_offset_v
      real(DP)                                :: mtot
      integer(I4B)                            :: i, nfrag
         
      mtot = sum(m_frag(:))
      nfrag = size(m_frag(:))
      mx_frag(:) = 0.0_DP
      mv_frag(:) = 0.0_DP

      do i = 1, nfrag
         mx_frag = mx_frag(:) + (x_frag(:,i) + xcom(:)) * m_frag(i)
         mv_frag = mv_frag(:) + (v_frag(:,i) + vcom(:)) * m_frag(i)
      end do
      COM_offset_x(:) = xcom(:) - mx_frag(:) / mtot
      COM_offset_v(:) = vcom(:) - mv_frag(:) / mtot
      do i = 1, nfrag 
         x_frag(:, i) = x_frag(:, i) + COM_offset_x(:)
         v_frag(:, i) = v_frag(:, i) + COM_offset_v(:)
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
      integer(I4B), intent(in) :: npl
      type(symba_pl), intent(in) :: symba_plA
      logical, dimension(:), intent(in) :: lexclude
      real(DP), intent(out) :: ke_orbit, ke_spin, pe
      real(DP), dimension(:), intent(out)      :: Ltot
      integer(I4B), intent(in), optional :: nfrag
      real(DP), dimension(:), intent(in), optional :: m_frag, rad_frag
      real(DP), dimension(:,:), intent(in), optional :: Ip_frag, xb_frag, vb_frag, rot_frag
      ! Internals
      integer(I4B) :: i, ntot, npl_loc, nplm
      type(symba_pl) :: symba_pl_loc
      real(DP) :: te

      npl_loc = count(.not.lexclude(:))
      ntot = npl_loc
      if (present(nfrag)) ntot = npl_loc + nfrag
      call symba_pl_allocate(symba_pl_loc, ntot)

      symba_pl_loc%helio%swiftest%mass(1:npl_loc) = pack(symba_plA%helio%swiftest%mass(1:npl), .not.lexclude(:))
      symba_pl_loc%helio%swiftest%radius(1:npl_loc) = pack(symba_plA%helio%swiftest%radius(1:npl), .not.lexclude(:))
      do i = 1, NDIM
         symba_pl_loc%helio%swiftest%xb(i,1:npl_loc) = pack(symba_plA%helio%swiftest%xb(i,1:npl), .not.lexclude(:))
         symba_pl_loc%helio%swiftest%vb(i,1:npl_loc) = pack(symba_plA%helio%swiftest%vb(i,1:npl), .not.lexclude(:))
         symba_pl_loc%helio%swiftest%rot(i,1:npl_loc) = pack(symba_plA%helio%swiftest%rot(i,1:npl), .not.lexclude(:))
         symba_pl_loc%helio%swiftest%Ip(i,1:npl_loc) = pack(symba_plA%helio%swiftest%Ip(i,1:npl), .not.lexclude(:))
      end do
      if (present(nfrag)) then
         symba_pl_loc%helio%swiftest%Ip(:,npl_loc+1:ntot) = Ip_frag(:,:)
         symba_pl_loc%helio%swiftest%mass(npl_loc+1:ntot) = m_frag(:)
         symba_pl_loc%helio%swiftest%radius(npl_loc+1:ntot) = m_frag(:)
         symba_pl_loc%helio%swiftest%xb(:,npl_loc+1:ntot) =  xb_frag(:,:)
         symba_pl_loc%helio%swiftest%vb(:,npl_loc+1:ntot) =  vb_frag(:,:)
         symba_pl_loc%helio%swiftest%rot(:,npl_loc+1:ntot) = rot_frag(:,:)
      end if

      symba_pl_loc%helio%swiftest%status(:) = ACTIVE
      
      call coord_b2h(ntot, symba_pl_loc%helio%swiftest)

      nplm = count(symba_pl_loc%helio%swiftest%mass > param%mtiny)
      call util_dist_index_plpl(ntot, nplm, symba_pl_loc)

      call symba_energy_eucl(ntot, symba_pl_loc, param%j2rp2, param%j4rp4, ke_orbit, ke_spin, pe, te, Ltot)

      call symba_pl_deallocate(symba_pl_loc)

      return
   end subroutine symba_frag_pos_energy_calc


end subroutine symba_frag_pos
