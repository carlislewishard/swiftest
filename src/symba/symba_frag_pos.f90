subroutine symba_frag_pos (symba_plA, idx_parents, x, v, L_spin, Ip, mass, radius, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
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

   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:), intent(in)        :: mass, radius, m_frag, rad_frag
   real(DP), dimension(:,:), intent(in)      :: Ip_frag
   real(DP), dimension(:,:), intent(out)     :: x_frag, v_frag, rot_frag
   integer(I4B), dimension(2), intent(inout) :: idx_parents
   type(symba_pl)                            :: symba_plA

   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, j, nfrag, fam_size, istart, non_fam_size, npl
   real(DP), dimension(NDIM)               :: v_cross_x, delta_v, delta_x, xcom, vcom
   real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm, KE_residual
   real(DP)                                :: f_anelastic, Etot_before, Etot_after, KE_before, U_before, v1mag2, v2mag2, U_p1_before, U_p2_before 
   real(DP)                                :: U_after, KE_spin_before, KE_spin_after, KE_after, KE_corrected, U_corrected, U_frag_after
   real(DP), dimension(NDIM)               :: v_col_unit_vec, tri_pro, tri_pro_unit_vec, v_com, vc1, vc2
   real(DP), dimension(NDIM)               :: Ltot, h
   real(DP)                                :: rot2, v2, ip_family, f_corrected, A, B
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9
   integer(I4B), dimension(:), allocatable :: family, non_family

   associate(nchild1 => symba_plA%kin(idx_parents(1))%nchild, nchild2 => symba_plA%kin(idx_parents(2))%nchild, &
             xhpl => symba_plA%helio%swiftest%xh, xbpl => symba_plA%helio%swiftest%xh, vbpl => symba_plA%helio%swiftest%vb, &
             Mpl => symba_plA%helio%swiftest%mass, Ippl => symba_plA%helio%swiftest%Ip, radpl => symba_plA%helio%swiftest%radius, &
             rotpl => symba_plA%helio%swiftest%rot, status => symba_plA%helio%swiftest%status)

      ! Find the center of mass of the collisional system
      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 

      ! Make the list of family members (bodies involved in the collision)
      fam_size = 2 + nchild1 + nchild2
      allocate(family(fam_size))
      family(1) = idx_parents(1)
      family(2) = idx_parents(2)
      istart = 2 + nchild1

      if (nchild1 > 0) family(3:istart) = symba_plA%kin(idx_parents(1))%child(1:nchild1)
      if (nchild2 > 0) family(istart+1:istart+1+nchild2) = symba_plA%kin(idx_parents(2))%child(1:nchild2)

      ! Make the list of non-family members (bodies not involved in the collision)
      npl = count(status(:) /= INACTIVE)
      non_fam_size = npl - fam_size
      allocate(non_family(non_fam_size))
      i = 0
      do j = 1, size(status(:))
         if (any(family(:) == j) .or. (status(j) == INACTIVE)) cycle
         i = i + 1
         non_family(i) = j
      end do

      U_before = 0.0_DP
      KE_before = 0.0_DP
      KE_spin_before = 0.0_DP
      Ltot = 0.0_DP

      r_col_norm = 0.0_DP
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
         Ltot(:) = Ltot(:) + Mpl(family(i)) * (h(:) + Ippl(3, family(i)) * radpl(family(i))**2 * rotpl(:,family(i)))
         KE_before = KE_before + 0.5_DP * Mpl(family(i)) * v2 
         KE_spin_before = KE_spin_before + 0.5_DP * Mpl(family(i)) * Ippl(3,family(i)) * rot2 * radpl(family(i))**2
         r_col_norm = r_col_norm + radpl(family(i)) !Mpl(family(i)) * norm2(xbpl(:,family(i)) - xcom(:))
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

      ! Calculate the triple product to get the plane of the fragment distribution
      delta_v(:) = v(:, 2) - v(:, 1)
      delta_x(:) = x(:, 2) - x(:, 1)
      call util_crossproduct(delta_v,delta_x,v_cross_x)
      call util_crossproduct(v_cross_x,delta_v,tri_pro)
      tri_pro_unit_vec(:) = tri_pro(:) / norm2(tri_pro(:))
      v_col_norm = norm2(delta_v(:))               ! pre-collision velocity magnitude
      v_col_unit_vec(:) = delta_v(:) / v_col_norm  ! unit vector of collision velocity 
      r_col_norm = r_col_norm / nfrag
      do i = 1, nfrag
         ! Place the fragments on the collision plane at a distance proportional to mass wrt the collisional barycenter
         ! This gets updated later after the new potential energy is calculated
         r_frag_norm = r_col_norm * mtot / m_frag(i)

         x_frag(:,i) =  r_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)  + &
                                       (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:)) 
                        
         ! Apply a simple mass weighting first to ensure that the velocity follows the barycenter
         ! This gets updated later after the new potential and kinetic energy is calcualted
         v_frag_norm = v_col_norm * mtot / m_frag(i)   
         v_frag(:,i) =  v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:) + &
                                       (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:)) 
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

      ! Adjust the fragment velocities so that they have the their total energy reduced by an amount set by the anelastic parameter
      ! Make sure we don't end up with negative energy (bound system). If so, we'll adjust the radius so that the potential energy
      ! takes up the negative part
      f_anelastic = 1.000_DP ! TODO: Should this be set by the user or kept as a constant?
      KE_residual = KE_after + KE_spin_after + U_after - f_anelastic * Etot_before  

      write(*,*) "SYMBA_FRAG_POS Etot_before : ", Etot_before
      write(*,*) "SYMBA_FRAG_POS Etot_after  : ", KE_after + KE_spin_after + U_after
      write(*,*) 'KE_residual: ',KE_residual

      A = 0.0_DP
      B = 0.0_DP

      do i = 1, nfrag
         A = A + m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
         B = B + m_frag(i) * dot_product(v_frag(:,i), vcom(:))
      end do

      if ((B + A)**2 - 2 * A * KE_residual > 0.0_DP) then
         f_corrected = (- B + sqrt((B + A)**2 - 2 * A * KE_residual)) / A
      else
         f_corrected = 0.0_DP
      end if
      !write(*,*) 'A: ',A
      !write(*,*) 'B: ',B
      write(*,*) 'f_corrected: ',f_corrected
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

      write(*,*) "SYMBA_FRAG_POS KE_corr   : ", KE_after + KE_spin_after
      write(*,*) "SYMBA_FRAG_POS KE_ratio  : ", (KE_after + KE_spin_after) / (KE_before + KE_spin_before)
      write(*,*) "SYMBA_FRAG_POS U_before  : ", U_before
      write(*,*) "SYMBA_FRAG_POS U_after   : ", U_after
      write(*,*) "SYMBA_FRAG_POS U_ratio   : ", U_after / U_before
      write(*,*) "SYMBA_FRAG_POS E_before  : ", KE_before + KE_spin_before + U_before
      write(*,*) "SYMBA_FRAG_POS E_after   : ", KE_after + KE_spin_after + U_after
      write(*,*) "SYMBA_FRAG_POS E_before / E_after   : ", Etot_before / Etot_after

      deallocate(family, non_family)
   end associate
   return 
end subroutine symba_frag_pos
