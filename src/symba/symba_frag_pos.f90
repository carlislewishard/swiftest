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
   real(DP), dimension(NDIM)               :: xc, vc, x_cross_v, delta_r, delta_v, xcom, vcom
   real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm, KE_residual
   real(DP)                                :: f_anelastic, Etot_before, Etot_after, KE_before, U_before
   real(DP)                                :: U_after, KE_spin_before, KE_spin_after, KE_after, KE_corrected, U_corrected, U_frag_after
   real(DP), dimension(NDIM)               :: r_col_unit_vec, v_col_unit_vec, v_plane_unit_vec
   real(DP), dimension(NDIM)               :: Ltot, h, dx
   real(DP)                                :: rot2, v2, f_corrected, A, B, C
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9
   integer(I4B), dimension(:), allocatable :: family, non_family
   real(DP)                                :: Esys, rmag


   associate(nchild1 => symba_plA%kin(idx_parents(1))%nchild, nchild2 => symba_plA%kin(idx_parents(2))%nchild, &
             xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status)

      npl = size(status(:))
      !****************************************************************k*
      ! Find the total system energy for reporting (testing only)
      Esys = 0.0_DP
      do i = 1, npl
         if (status(i) == INACTIVE) cycle
         v2 = dot_product(vbpl(:,i), vbpl(:,i))
         rot2 = dot_product(rotpl(:,i), rotpl(:,i))
         Esys = Esys + 0.5_DP * Mpl(i) * (v2 + Ippl(3,i) * radpl(i)**2 * rot2)
      end do
      do i = 1, npl - 1
         if (status(i) == INACTIVE) cycle
         do j = i + 1, npl
            if (status(j) == INACTIVE) cycle
            dx(:) = xbpl(:, j) - xbpl(:, i) 
            rmag = norm2(dx(:)) 
            if (rmag > tiny(rmag)) Esys = Esys - Mpl(i) * Mpl(j) / rmag 
         end do
      end do
      !****************************************************************k*

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
      family(1) = idx_parents(1)
      family(2) = idx_parents(2)
      istart = 2 + nchild1

      if (nchild1 > 0) family(3:istart) = symba_plA%kin(idx_parents(1))%child(1:nchild1)
      if (nchild2 > 0) family(istart+1:istart+1+nchild2) = symba_plA%kin(idx_parents(2))%child(1:nchild2)
      fam_size = count(status(family(:)) == ACTIVE)
      family(:) = pack(family(:), status(family(:)) == ACTIVE)

      ! Make the list of non-family members (bodies not involved in the collision)
      non_fam_size = count(status(:) /= INACTIVE) - fam_size
      allocate(non_family(non_fam_size))
      i = 0
      do j = 1, npl
         if (any(family(:) == j) .or. (status(j) == INACTIVE)) cycle
         i = i + 1
         non_family(i) = j
      end do

      U_before = 0.0_DP
      KE_before = 0.0_DP
      KE_spin_before = 0.0_DP

      do i = 1, fam_size
         ! Calculate the potential energy between family members
         do j = i + 1, fam_size
            U_before = U_before - Mpl(family(i)) * Mpl(family(j)) / norm2(xbpl(:,family(i)) - xbpl(:,family(j)))
         end do
         
         ! Add the contribution due to non-family members
         do j = 1, non_fam_size
            U_before = U_before - Mpl(family(i)) * Mpl(non_family(j)) / norm2(xbpl(:,family(i)) - xbpl(:,non_family(j)))
         end do

         !! Calculate the orbital and spin kinetic energy
         v2   = dot_product(vbpl(:,family(i)), vbpl(:,family(i)))
         rot2 = dot_product(rotpl(:,family(i)), rotpl(:,family(i)))
         call util_crossproduct(xbpl(:,family(i)), vbpl(:,family(i)), h(:))
         
         !! Calculate the orbital and rotational kinetic energy between the two bodies 
         KE_before = KE_before + 0.5_DP * Mpl(family(i)) * v2 
         KE_spin_before = KE_spin_before + 0.5_DP * Mpl(family(i)) * Ippl(3,family(i)) * rot2 * radpl(family(i))**2
      end do

      Etot_before = KE_before + KE_spin_before + U_before

      ! Now create the fragment distribution
      nfrag = size(x_frag, 2)

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

      U_after = 0.0_DP
      ! Calculate the new potential energy of the system 
      do i = 1, nfrag
         ! Calculate the contribution of the potential energy between fragments
         do j = i+1,nfrag
            U_after = U_after - m_frag(i) * m_frag(j) / norm2(x_frag(:,i) - x_frag(:,j))
         end do
         ! Add the contribution due to non-family members
         do j = 1, non_fam_size
            U_after = U_after - m_frag(i) * Mpl(non_family(j)) / norm2(x_frag(:,i) + xcom(:) - xbpl(:,non_family(j)))
         end do
      end do

      ! Compute the current kinetic energy of the system so we can scale the velocity vectors to the correct value
      KE_after = 0.0_DP
      KE_spin_after = 0.0_DP
      do i = 1, nfrag
         KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i) + vcom(:), v_frag(:,i) + vcom(:))
         KE_spin_after = KE_spin_after + 0.5_DP * m_frag(i) * rad_frag(i)**2 * Ip_frag(3,i) * dot_product(rot_frag(:,i), rot_frag(:,i))
      end do
      Etot_after = KE_after + KE_spin_after + U_after

      ! Adjust the fragment velocities so that they have the their total energy reduced by an amount set by the anelastic parameter
      ! Make sure we don't end up with negative energy (bound system). If so, we'll adjust the radius so that the potential energy
      ! takes up the negative part

      !f_anelastic = 0.9999999_DP ! TODO: Should this be set by the user or kept as a constant?
      !KE_residual = KE_spin_after + U_after - U_before - f_anelastic * (KE_before + KE_spin_before) 

      KE_residual = KE_spin_after + U_after - U_before - KE_before - KE_spin_before + Qloss

100 format (A14,5(ES9.2,1X,:))
      write(*,   "('              Energy normalized by |Esystem_original|')")
      write(*,   "('             |    T_orb    T_spin         T         U      Etot')")
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' original    |',KE_before / abs(Esys), KE_spin_before / abs(Esys),(KE_before + KE_spin_before) / abs(Esys), U_before / abs(Esys),Etot_before/abs(Esys)
      write(*,100) ' first pass  |',KE_after / abs(Esys), KE_spin_after / abs(Esys), (KE_after + KE_spin_after) / abs(Esys), U_after / abs(Esys), Etot_after / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' change      |',(KE_after - KE_before) / abs(Esys), (KE_spin_after - KE_spin_before)/ abs(Esys), (KE_after + KE_spin_after - KE_before - KE_spin_before)/ abs(Esys), (U_after - U_before) / abs(Esys), (Etot_after-Etot_before) / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' T_res       |',KE_residual / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' Q_loss      |',Qloss / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")

      A = 0.0_DP
      B = 0.0_DP
      C = 2 * KE_residual

      do i = 1, nfrag
         A = A + m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
         B = B + m_frag(i) * dot_product(v_frag(:,i), vcom(:))
         C = C + m_frag(i) * dot_product(vcom(:), vcom(:))
      end do

      if ((B**2 - A * C) > 0.0_DP) then
         f_corrected = (- B + sqrt(B**2 - A * C)) / A
         lmerge = .false.
      else
         f_corrected = 0.0_DP
         lmerge = .true.
      end if
      v_frag(:,:) =  f_corrected * v_frag(:,:) 

      ! Shift the fragments into the system barycenter frame
      do i = 1, nfrag
         x_frag(:,i) = x_frag(:, i) + xcom(:)
         v_frag(:,i) = v_frag(:, i) + vcom(:)
      end do

      ! REMOVE THE FOLLOWING AFTER TESTING
      ! Calculate the new energy of the system of fragments
      KE_after = 0.0_DP
      do i = 1, nfrag
         KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
      end do
      Etot_after = KE_after + KE_spin_after + U_after
      write(*,100) ' f_corrected |',f_corrected
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' final       |',KE_after / abs(Esys), KE_spin_after / abs(Esys), (KE_after + KE_spin_after) / abs(Esys), U_after / abs(Esys), Etot_after / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,100) ' change      |',(KE_after - KE_before) / abs(Esys), (KE_spin_after - KE_spin_before)/ abs(Esys), (KE_after + KE_spin_after - KE_before - KE_spin_before)/ abs(Esys), (U_after - U_before) / abs(Esys), (Etot_after-Etot_before) / abs(Esys)
      write(*,   "(' ------------------------------------------------------------------')")
      write(*,*)   

      deallocate(family, non_family)
   end associate
   return 
end subroutine symba_frag_pos
