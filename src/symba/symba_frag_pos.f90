subroutine symba_frag_pos (symba_plA, idx_parents, x, v, L_spin, Ip, mass, radius, &
                           Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag, lmerge, Qloss)
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

   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(:), intent(in)    :: idx_parents
   real(DP), intent(in)                      :: Qloss
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:), intent(in)        :: mass, radius, m_frag, rad_frag
   real(DP), dimension(:,:), intent(in)      :: Ip_frag
   real(DP), dimension(:,:), intent(out)     :: x_frag, v_frag, rot_frag
   logical, intent(out)                      :: lmerge ! Answers the question: Should this have been a merger instead?

   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, j, k, nfrag, fam_size, istart, non_fam_size, npl
   real(DP), dimension(NDIM)               :: xc, vc, x_cross_v, delta_r, delta_v, xcom, vcom, v_phi_unit, v_r_unit
   real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm, KE_residual
   real(DP)                                :: f_anelastic, Etot_before, Etot_after, KE_before, U_before
   real(DP)                                :: U_after, KE_spin_before, KE_spin_after, KE_after, KE_corrected, U_corrected, U_frag_after, KE_family
   real(DP), dimension(NDIM)               :: r_col_unit_vec, v_col_unit_vec, v_plane_unit_vec
   real(DP), dimension(NDIM)               :: Ltot, h, dx
   real(DP)                                :: rot2, v2, f_corrected, A, B, C
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9
   integer(I4B), dimension(:), allocatable :: family, non_family
   real(DP)                                :: rmag
   logical, dimension(:), allocatable      :: lfamily
   real(DP), dimension(:,:), allocatable   :: xtmp, vtmp, v_phi, v_r
   logical, dimension(:), allocatable      :: lexclude
   
   associate(nchild1 => symba_plA%kin(idx_parents(1))%nchild, nchild2 => symba_plA%kin(idx_parents(2))%nchild, &
             xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status, npl => symba_plA%helio%swiftest%nbody)

      allocate(lexclude(npl))
      lexclude(:) = status(1:npl) == INACTIVE
      call symba_frag_pos_energy(npl, symba_plA, lexclude, KE_before, KE_spin_before, U_before)
      Etot_before = KE_before + KE_spin_before + U_before

      ! Find the center of mass of the collisional system
      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 
      
      ! Relative position and velocity vectors of the two impacting "clouds" 
      delta_r(:) = x(:, 2) - x(:, 1)
      r_col_norm = norm2(delta_r(:))
      delta_v(:) = v(:, 2) - v(:, 1)
      v_col_norm = norm2(delta_v(:))               

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
      KE_family = 0.0_DP
      do i = 1, fam_size
         KE_family = KE_family + Mpl(family(i)) * dot_product(vbpl(:,family(i)), vbpl(:,family(i)))
      end do
      KE_family = 0.5_DP * KE_family

      ! Make the list of non-family members (bodies not involved in the collision)
      non_fam_size = count(status(:) /= INACTIVE) - fam_size
      allocate(non_family(non_fam_size))
      i = 0
      do j = 1, npl
         if (any(family(:) == j)) lexclude(j) = .true.
         if (any(family(:) == j) .or. (status(j) == INACTIVE)) cycle
         i = i + 1
         non_family(i) = j
      end do

      ! Now create the fragment distribution
      nfrag = size(x_frag, 2)
      allocate(v_r(NDIM,nfrag))
      allocate(v_phi(NDIM,nfrag))

      ! Calculate the position of each fragment 
      ! Theta is a phase shift value that ensures that successive nearby collisions in a single step are rotated to avoid possible overlap
      theta = (2 * PI) / nfrag
      ! Shifts the starting circle of fragments around so that multiple fragments generated 
      ! from a single collision in a single time step don't pile up on top of each other
      phase_ang = theta * thetashift / SHIFTMAX
      thetashift = thetashift + 1
      IF (thetashift >= shiftmax) thetashift = 0

      ! Compute orbital angular momentum of pre-impact system. This will be the normal vector to the collision fragment plane
      Ltot = L_spin(:,1) + L_spin(:,2)
      do j = 1, 2
         xc(:) = x(:, j) - xcom(:)
         vc(:) = v(:, j) - vcom(:)
         call utiL_crossproduct(xc(:), vc(:), x_cross_v(:))
         Ltot(:) = Ltot(:) + mass(j) * x_cross_v(:)
      end do

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
      ! Adjust fragment positionx so that their center of mass follows the original pre-impact trajectory
      call symba_frag_adjust (xcom, vcom, x, v, mass, radius, L_spin, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)

      allocate(xtmp, source=x_frag)
      allocate(vtmp, source=v_frag)
      do i = 1, nfrag
         xtmp(:,i) = xtmp(:,i) + xcom(:)
         vtmp(:,i) = vtmp(:,i) + vcom(:)
      end do
      call symba_frag_pos_energy(npl, symba_plA, lexclude, KE_after, KE_spin_after, U_after, &
         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, x_frag=xtmp, v_frag=vtmp, rot_frag=rot_frag)
      Etot_after = KE_after + KE_spin_after + U_after

      ! Adjust the fragment velocities so that they have the their total energy reduced by an amount set by the anelastic parameter
      ! Make sure we don't end up with negative energy (bound system). If so, we'll adjust the radius so that the potential energy
      ! takes up the negative part

      KE_residual = KE_family + (KE_spin_before - KE_spin_after) + (U_before - U_after) - Qloss

100 format (A14,5(ES9.2,1X,:))
      write(*,   "('              Energy normalized by |Etot_before|')")
      write(*,   "('             |    T_orb    T_spin         T         U      Etot')")
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' original    |',KE_before / abs(Etot_before), KE_spin_before / abs(Etot_before),(KE_before + KE_spin_before) / abs(Etot_before), U_before / abs(Etot_before),Etot_before/abs(Etot_before)
      write(*,100) ' first pass  |',KE_after / abs(Etot_before), KE_spin_after / abs(Etot_before), (KE_after + KE_spin_after) / abs(Etot_before), U_after / abs(Etot_before), Etot_after / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' change      |',(KE_after - KE_before) / abs(Etot_before), (KE_spin_after - KE_spin_before)/ abs(Etot_before), (KE_after + KE_spin_after - KE_before - KE_spin_before)/ abs(Etot_before), (U_after - U_before) / abs(Etot_before), (Etot_after-Etot_before) / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' T_res       |',KE_residual / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' Q_loss      |',Qloss / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")

      do i = 1, nfrag
         v_r(:,i) = -dot_product(v_frag(:,i), x_frag(:,i)) * x_frag(:,i) / dot_product(x_frag(:,i), x_frag(:,i))
         v_phi(:,i) = v_frag(:,i) - v_r(:,i) 
      end do

      A = 0.0_DP
      B = 0.0_DP

      do i = 1, nfrag
         A = A + m_frag(i) * dot_product(v_r(:,i), v_r(:,i))
         C = m_frag(i) * (dot_product(vcom(:), vcom(:)) + dot_product(v_phi(:,i), v_phi(:,i)))
      end do
      C = 2 * KE_residual - C

      B = C / A
      if (B > 0.0_DP) then
         f_corrected = sqrt(B)
         lmerge = .false.
      else
         f_corrected = 0.0_DP
         lmerge = .true.
      end if

      ! Shift the fragments into the system barycenter frame
      do i = 1, nfrag
         x_frag(:,i) = x_frag(:, i) + xcom(:)
         v_frag(:,i) = f_corrected * v_r(:, i) + v_phi(:, i) + vcom(:)
      end do

      ! REMOVE THE FOLLOWING AFTER TESTING
      ! Calculate the new energy of the system of fragments
      call symba_frag_pos_energy(npl, symba_plA, lexclude, KE_after, KE_spin_after, U_after, &
         nfrag=nfrag, Ip_frag=Ip_frag, m_frag=m_frag, rad_frag=rad_frag, x_frag=x_frag, v_frag=v_frag, rot_frag=rot_frag)
      Etot_after = KE_after + KE_spin_after + U_after
      write(*,100) ' f_corrected |',f_corrected
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' final       |',KE_after / abs(Etot_before), KE_spin_after / abs(Etot_before), (KE_after + KE_spin_after) / abs(Etot_before), U_after / abs(Etot_before), Etot_after / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' change      |',(KE_after - KE_before) / abs(Etot_before), (KE_spin_after - KE_spin_before)/ abs(Etot_before), (KE_after + KE_spin_after - KE_before - KE_spin_before)/ abs(Etot_before), (U_after - U_before) / abs(Etot_before), (Etot_after-Etot_before) / abs(Etot_before)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,*)   

      deallocate(family, non_family)
   end associate
   return 

   contains

   subroutine symba_frag_pos_energy(npl, symba_plA, lexclude, KE, KE_spin, U, nfrag, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
      implicit none
      integer(I4B), intent(in) :: npl
      type(symba_pl), intent(in) :: symba_plA
      logical, dimension(:), intent(in) :: lexclude
      real(DP), intent(out) :: KE, KE_spin, U
      integer(I4B), intent(in), optional :: nfrag
      real(DP), dimension(:), intent(in), optional :: m_frag, rad_frag
      real(DP), dimension(:,:), intent(in), optional :: Ip_frag, x_frag, v_frag, rot_frag
      integer(I4B) :: i
      real(DP), dimension(NDIM) :: dx
      real(DP) :: rmag

      associate(xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
         Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
         rotpl => symba_plA%helio%swiftest%rot)
         U = 0.0_DP
         KE = 0.0_DP
         KE_spin = 0.0_DP
         do i = 1, npl
            if (lexclude(i)) cycle
            v2 = dot_product(vbpl(:,i), vbpl(:,i))
            rot2 = dot_product(rotpl(:,i), rotpl(:,i))
            KE = KE + Mpl(i) * v2 
            KE_spin = KE_spin + Mpl(i) * Ippl(3,i) * radpl(i)**2 * rot2
         end do

         do i = 1, npl - 1
            if (lexclude(i)) cycle
            do j = i + 1, npl
               if (lexclude(j)) cycle
               dx(:) = xbpl(:, j) - xbpl(:, i) 
               rmag = norm2(dx(:)) 
               if (rmag > tiny(rmag)) U = U - Mpl(i) * Mpl(j) / rmag 
            end do
            if (present(nfrag)) then
               do j = 1, nfrag
                  dx(:) = x_frag(:, j) - xbpl(:, i) 
                  rmag = norm2(dx(:)) 
                  if (rmag > tiny(rmag)) U = U - Mpl(i) * m_frag(j) / rmag 
               end do
            end if
         end do
         if (present(nfrag)) then
            do i = 1, nfrag
               v2 = dot_product(v_frag(:,i), v_frag(:,i))
               rot2 = dot_product(rot_frag(:,i), rot_frag(:,i))
               KE = KE + m_frag(i) * v2 
               KE_spin = KE_spin + m_frag(i) * Ip_frag(3,i) * rad_frag(i)**2 * rot2
               do j = i + 1, nfrag
                  dx(:) = x_frag(:, j) - x_frag(:, i) 
                  rmag = norm2(dx(:)) 
                  if (rmag > tiny(rmag)) U = U - m_frag(i) * m_frag(j) / rmag
               end do
            end do
         end if

         KE = 0.5_DP * KE
         KE_spin = 0.5_DP * KE_spin
      end associate

      return
   end subroutine symba_frag_pos_energy
end subroutine symba_frag_pos
