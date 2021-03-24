subroutine symba_frag_pos (x, v, L_spin, Ip, mass, radius, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
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

   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, j, nfrag
   real(DP), dimension(NDIM)               :: v_cross_x, delta_v, delta_x, L_spin_tot
   real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
   real(DP)                                :: f_anelastic, Etot_before, KE_before, U_before, v1mag2, v2mag2
   real(DP)                                :: U_after, KE_spin_after, KE_after
   real(DP), dimension(NDIM)               :: v_col_unit_vec, tri_pro, tri_pro_unit_vec, v_com, vc1, vc2
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9

   nfrag = size(x_frag, 2)
   mtot = sum(mass(:))

   ! Now work out where to put the new bodies such that we conserve momentum and lose an appropriate amount of energy
   ! Find collision velocity
   delta_v(:) = v(:, 2) - v(:, 1)
   delta_x(:) = x(:, 2) - x(:, 1)

   v_com(:) = ((v(:, 1) * mass(1)) + (v(:, 2) * mass(2))) / mtot

   v_col_norm = norm2(delta_v(:)) ! pre-collision velocity magnitude
   r_col_norm = norm2(delta_x(:)) ! pre-collision distance 
   v_col_unit_vec(:) = delta_v(:) / v_col_norm ! unit vector of collision velocity 

   ! Determine the radius and angular spacing of fragments placed in a circle around the COM of the collision
   !r_circle = sum(rhill(:)) / (2 * sin(PI / nfrag))
   theta = (2 * PI) / nfrag
   ! Shifts the starting circle of fragments around so that multiple fragments generated 
   ! from a single collision in a single time step don't pile up on top of each other
   phase_ang = theta * thetashift / SHIFTMAX
   thetashift = thetashift + 1
   IF (thetashift >= shiftmax) thetashift = 0

   ! Fragment velocity is the  collision velocity with an adjustment for new distance using vis viva
   !v_frag_norm =  sqrt(v_col_norm**2 - 2 * mtot * (1._DP / r_col_norm - 1._DP / r_circle))

   ! Calculate the triple product to get the plane of the fragment distribution
   call util_crossproduct(delta_v,delta_x,v_cross_x)
   call util_crossproduct(v_cross_x,delta_v,tri_pro)

   tri_pro_unit_vec(:) = tri_pro(:) / norm2(tri_pro(:))

   f_anelastic = 0.1_DP ! TODO: Should this be set by the user or kept as a constant?
   vc1(:) = v(:,1) - v_com(:)
   vc2(:) = v(:,2) - v_com(:)
   v1mag2 = dot_product(vc1(:), vc1(:))
   v2mag2 = dot_product(vc2(:), vc2(:))

   KE_before = 0.5_DP * (mass(1) * v1mag2 + mass(2) * v2mag2) 

   ! Add in kinetic energy from spin
   rot(:,1) = L_spin(:,1) / (mass(1) * Ip(3,1) * radius(1)**2)
   rot(:,2) = L_spin(:,2) / (mass(2) * Ip(3,2) * radius(2)**2)
   L_spin_tot(:) = L_spin(:,1) + L_spin(:,2)
   KE_before = KE_before + 0.5_DP * (dot_product(rot(:,1),L_spin(:,1)) + dot_product(rot(:,2),L_spin(:,2)))

   ! Add in potential energy
   U_before = -mass(1) * mass(2) / r_col_norm
   Etot_before = KE_before + U_before
   ! Calculate the position of each fragment 
   KE_spin_after = 0.0_DP
   do i = 1, nfrag
      ! Place the fragments on the collision plane at a distance proportional to mass wrt the collisional barycenter
      r_frag_norm = r_col_norm * mtot / m_frag(i) 
      x_frag(:,i) =  r_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)  + &
                                    (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))
      ! Assume spin angular momentum is conserved to get the rotation rates of the fragments
      ! This will get updated later in symba_frag_adjust to conserve total system angular momentum
      rot_frag(:,i) = L_spin_tot(:) / nfrag / (m_frag(i) * rad_frag(i)**2 * Ip_frag(3,i))
      KE_spin_after = KE_spin_after + 0.5_DP * dot_product(rot_frag(:,i), L_spin_tot(:) / nfrag)
   end do

   ! Calculate the new potential energy of the system of fragments
   U_after = 0.0_DP
   do j = 1, nfrag - 1
      do i = j+1, nfrag
         delta_x(:) = x_frag(:,i) - x_frag(:,j)
         U_after = U_after - m_frag(i) * m_frag(j) / norm2(delta_x(:))
      end do 
   end do

   ! Adjust fragment positions so that they have the same potential energy as the original collisional pair
   x_frag(:,:) = x_frag(:,:) * U_after / U_before

   KE_after = KE_spin_after
   do i=1, nfrag 
      ! Include the spin kinetic energy when computing the new kinetic energy of the fragment
      v_frag_norm = sqrt(2 * (f_anelastic * Etot_before - KE_spin_after - U_before) / (nfrag * m_frag(i))) 
      v_frag(:,i) =  v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:) + &
                                 (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))
      KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
   end do


   return 
end subroutine symba_frag_pos
