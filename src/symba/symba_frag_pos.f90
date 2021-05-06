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
   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, j, nfrag, fam_size, istart
   real(DP), dimension(NDIM)               :: xcom, vcom, Ltot
   real(DP)                                :: mtot, Ltot_before, Ltot_after
   real(DP)                                :: Etot_before, Etot_after, ke_before, pe_before
   real(DP)                                :: pe_after, ke_spin_before, ke_spin_after, ke_after, ke_family, ke_target, ke_frag
   real(DP), dimension(NDIM)               :: h, dx
   real(DP)                                :: rmag
   logical, dimension(:), allocatable      :: lexclude
   character(len=*), parameter             :: fmtlabel = "(A14,10(ES9.2,1X,:))"
   real(DP), dimension(NDIM)               :: x_cross_v, v_phi_unit, h_unit, v_r_unit
   real(DP), dimension(:,:), allocatable   :: v_r, v_phi
   
   associate(xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status, npl => symba_plA%helio%swiftest%nbody, name => symba_plA%helio%swiftest%id)

      allocate(x_frag, source=xb_frag)
      allocate(v_frag, source=vb_frag)
      allocate(v_r, mold=v_frag)
      allocate(v_phi, mold=v_frag)
      fam_size = size(family)

      ! Find the center of mass of the collisional system
      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 

      allocate(lexclude(npl))
      where (status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
         lexclude(1:npl) = .true.  
      elsewhere
         lexclude(1:npl) = .false. 
      end where

      call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_before, ke_spin_before, pe_before, Ltot)
      Ltot_before = norm2(Ltot(:))
      Etot_before = ke_before + ke_spin_before + pe_before

      ! We need the original kinetic energy of just the pre-impact family members in order to balance the energy later
      ke_family = 0.0_DP
      do i = 1, fam_size
         ke_family = ke_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i))) !
         lexclude(family(i)) = .true. ! For all subsequent energy calculations the pre-impact family members will be replaced by the fragments
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

      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,        "('              Energy normalized by |Etot_before|')")
      !write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,        "('  First pass to get angular momentum ')")
      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
      !                                   (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
      !                                   (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
      !                                   (pe_after - pe_before) / abs(Etot_before), &
      !                                   (Etot_after - Etot_before) / abs(Etot_before), &
      !                                   (Ltot_after - Ltot_before) / Ltot_before
      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,        "('  Second pass to get energy ')")
      !write(*,        "(' ---------------------------------------------------------------------------')")
      !write(*,fmtlabel) ' Q_loss      |',-Qloss / abs(Etot_before)
      !write(*,        "(' ---------------------------------------------------------------------------')")

      
         ! Set the "target" ke_after (the value of the orbital kinetic energy that the fragments ought to have)
         ke_target = ke_family + (ke_spin_before - ke_spin_after) + (pe_before - pe_after) - Qloss
         call symba_frag_pos_kinetic_energy(xcom, vcom, m_frag, x_frag, v_frag, ke_target, lmerge)
         !write(*,        "(' ---------------------------------------------------------------------------')")
         !write(*,fmtlabel) ' T_frag targ |',ke_target / abs(Etot_before)

         ! Shift the fragments into the system barycenter frame
         do i = 1, nfrag
            xb_frag(:,i) = x_frag(:, i) + xcom(:)
            vb_frag(:,i) = v_frag(:, i) + vcom(:)
         end do

         ke_frag = 0.0_DP
         do i = 1, nfrag
            ke_frag = ke_frag + m_frag(i) * dot_product(vb_frag(:,i), vb_frag(:,i)) 
         end do
         ke_frag = 0.5_DP * ke_frag
         !write(*,fmtlabel) ' T_frag new  |',ke_frag / abs(Etot_before)
         !write(*,        "(' ---------------------------------------------------------------------------')")

         ! Calculate the new energy of the system of fragments
         call symba_frag_pos_energy_calc(npl, symba_plA, lexclude, ke_after, ke_spin_after, pe_after, Ltot,&
               nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, xb_frag=xb_frag, vb_frag=vb_frag, rot_frag=rot_frag)
         Etot_after = ke_after + ke_spin_after + pe_after
         Ltot_after = norm2(Ltot(:))
      
         lmerge = lmerge .or. ((Etot_after - Etot_before) / abs(Etot_before) > 0._DP) 

!      write(*,        "(' ---------------------------------------------------------------------------')")
!      write(*,        "('             |    T_orb    T_spin         T         pe      Etot      Ltot')")
!      write(*,        "(' ---------------------------------------------------------------------------')")
!      write(*,fmtlabel) ' change      |',(ke_after - ke_before) / abs(Etot_before), &
!                                       (ke_spin_after - ke_spin_before)/ abs(Etot_before), &
!                                       (ke_after + ke_spin_after - ke_before - ke_spin_before)/ abs(Etot_before), &
!                                       (pe_after - pe_before) / abs(Etot_before), &
!                                       (Etot_after - Etot_before) / abs(Etot_before), &
!                                       (Ltot_after - Ltot_before) / Ltot_before
!      write(*,        "(' ---------------------------------------------------------------------------')")
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
      real(DP)                                :: mtot, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
      real(DP)                                :: ecc_ellipse, b2a,  phase_ang, imp_param
      real(DP), dimension(NDIM)               :: Ltot, xc, vc, x_cross_v, delta_r, delta_v
      real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
      integer(I4B)                            :: i, nfrag
      real(DP), dimension(2,2)                :: orientation
      real(DP), dimension(2)                  :: frag_vec


      ! Now create the fragment distribution
      nfrag = size(m_frag(:))
      mtot = sum(mass(:))


      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      Ltot = L_spin(:,1) + L_spin(:,2)
      do i = 1, 2
         xc(:) = x(:, i) - xcom(:)
         vc(:) = v(:, i) - vcom(:)
         call utiL_crossproduct(xc(:), vc(:), x_cross_v(:))
         Ltot(:) = Ltot(:) + mass(i) * x_cross_v(:)
      end do
      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))     
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))

      ! We will initialize fragments on a planet defined by the pre-impact system, with the z-axis aligned with the angular momentum
      ! and the x-axis aligned with the impact velocity vector.
      z_col_unit = Ltot(:) / norm2(Ltot(:))
      x_col_unit(:) = delta_v(:) / v_col_norm  
      ! The cross product of the z- by x-axis will give us the y-axis
      call util_crossproduct(z_col_unit, x_col_unit, y_col_unit)


      
      ! Place the fragments on the collision planea on an ellipse, but with the distance proportional to mass wrt the collisional barycenter
      ! This gets updated later after the new potential energy is calculated
      ecc_ellipse = 0.90_DP


      b2a = 1.0_DP / sqrt(1.0_DP - ecc_ellipse**2)
      
      ! The orientation and angular spacing of fragments on the ellipse
      theta = (2 * PI) / nfrag
      ! Impirically determined phase angle that depends on the impact paarameter
      imp_param = norm2(Ltot(:)) / (r_col_norm * v_col_norm * mtot)
      phase_ang = 0.5_DP * PI 
      orientation = reshape([cos(phase_ang), sin(phase_ang), -sin(phase_ang), cos(phase_ang)], shape(orientation))

      ! Re-normalize position and velocity vectors by the fragment number so that for our initial guess we weight each
      ! fragment position by the mass and assume equipartition of energy for the velocity
      v_col_norm = 0.0_DP
      v_col_norm = v_col_norm / sqrt(1.0_DP * nfrag)
      r_col_norm = 2 * max(r_col_norm, sum(radius(:))) / nfrag ! To ensure that the new fragments aren't overlapping we will pick an initial starting radius 
      do i = 1, nfrag
         frag_vec(:) = [b2a * cos(theta * (i - 1)), sin(theta * (i - 1))]
         frag_vec(:) = matmul(orientation(:,:), frag_vec(:))

         r_frag_norm = r_col_norm * mtot / m_frag(i) 
         x_frag(:,i) =  r_frag_norm * (frag_vec(1) * x_col_unit(:) + frag_vec(2) * y_col_unit(:))
                        
         ! Apply a simple mass weighting first to ensure that the velocity follows the barycenter
         ! This gets updated later after the new potential and kinetic energy is calcualted
         v_frag_norm = v_col_norm * sqrt(mtot / m_frag(i))
         v_frag(:,i) = v_frag_norm * (frag_vec(1) * x_col_unit(:) + frag_vec(2) * y_col_unit(:))

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
      real(DP)                                :: f_corrected, mtot, A, C, rterm
      integer(I4B)                            :: i, nfrag
      real(DP), dimension(:,:), allocatable   :: v_r, v_phi
      real(DP), dimension(NDIM)               :: x_cross_v, v_phi_unit, h_unit, v_r_unit
         
      nfrag = size(m_frag)
      mtot = sum(m_frag)

      allocate(v_r(NDIM,nfrag))
      allocate(v_phi(NDIM,nfrag))

      do i = 1, nfrag
         call utiL_crossproduct(x_frag(:,i), v_frag(:,i), x_cross_v(:))
         h_unit(:) = x_cross_v(:) / norm2(x_cross_v(:))
         v_r_unit(:) = x_frag(:,i) / norm2(x_frag(:, i))
         call utiL_crossproduct(h_unit(:), v_r_unit(:), v_phi_unit(:))
         v_r(:,i) = v_r_unit(:) !dot_product(v_frag(:,i), v_r_unit(:)) * v_r_unit(:)
         v_phi(:,i) = dot_product(v_frag(:,i), v_phi_unit(:)) * v_phi_unit(:)
      end do

      C = 2 * ke_target
      do i = 1, nfrag
         C = C - m_frag(i) * (dot_product(vcom(:),vcom(:)) + dot_product(v_phi(:,i),v_phi(:,i)) + dot_product(vcom(:), v_phi(:, i)))
      end do

      if (C > 0.0_DP) then
         f_corrected = sqrt(C / mtot)
         lmerge = .false.
      else
         f_corrected = 0.0_DP
         lmerge = .true.
      end if

      ! Shift the fragments into the system barycenter frame
      do i = 1, nfrag
         v_frag(:,i) = f_corrected * v_r(:, i) + v_phi(:, i)
      end do

      call symba_frag_pos_com_adjust(xcom, vcom, m_frag, x_frag, v_frag)

      !write(*,fmtlabel) ' f_corrected |',f_corrected

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
         symba_plwksp%helio%swiftest%radius(npl+1:npl_new) = m_frag(:)
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
