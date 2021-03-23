subroutine symba_frag_pos (mtot, m1, m2, rhill, x, v, m_frag, x_frag, v_frag)
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

   real(DP), intent(in)                      :: mtot, m1, m2
   real(DP), dimension(:), intent(in)        :: rhill, m_frag
   real(DP), dimension(:,:), intent(in)      :: x, v
   real(DP), dimension(:,:), intent(out)     :: x_frag, v_frag

   integer(I4B)                            :: i, nfrag
   real(DP), dimension(NDIM)               :: v_cross_x, delta_v, delta_x
   real(DP)                                :: phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
   real(DP)                                :: f_anelastic, Etot_before
   real(DP), dimension(NDIM)               :: KE_before, U_before, v_com
   real(DP), dimension(NDIM)               :: v_col_unit_vec, tri_pro, tri_pro_unit_vec
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9

   nfrag = size(x_frag, 2)
   ! Now work out where to put the new bodies such that we conserve momentum and lose an appropriate amount of energy
   ! Find collision velocity
   delta_v(:) = v(:, 2) - v(:, 1)
   delta_x(:) = x(:, 2) - x(:, 1)

   v_com = ((v(:, 1) * m1) + (v(:, 2) * m2)) / (m1 + m2)

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

   ! Calculate the position and velocity of each fragment 
   do i=1, nfrag 
      f_anelastic = 0.1_DP
      KE_before(:) = (0.5_DP * m1 * (v(:,1) - v_com)**2) + (0.5_DP * m2 * (v(:,2) - v_com)**2)
      U_before(:) = (m1 * m2) / delta_x(:) 
      Etot_before = norm2(KE_before(:) + U_before(:))
      v_frag_norm = sqrt(((2.0_DP * f_anelastic) / (nfrag * m_frag(i))) * Etot_before)
      r_frag_norm = r_col_norm * ((m_frag(i) * mtot) / (m1 * m2))

      v_frag(:,i) =  v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:) + &
                                 (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))
      x_frag(:,i) =  r_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)  + &
                                 (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))

   end do
   return 
end subroutine symba_frag_pos
