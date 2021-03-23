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

   integer(I4B)                            :: i, j, nfrag
   real(DP), dimension(NDIM)               :: v_cross_x, delta_v, delta_x
   real(DP)                                :: phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
   real(DP)                                :: f_anelastic, Etot_before, KE_before, U_before, U_after, KE_after, v1mag2, v2mag2
   real(DP), dimension(NDIM)               :: v_col_unit_vec, tri_pro, tri_pro_unit_vec, v_com, vc1, vc2
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

   f_anelastic = 0.1_DP
   vc1(:) = v(:,1) - v_com(:)
   vc2(:) = v(:,2) - v_com(:)
   v1mag2 = dot_product(vc1(:), vc1(:))
   v2mag2 = dot_product(vc2(:), vc2(:))

   KE_before = 0.5_DP * (m1 * v1mag2 + m2 * v2mag2)
   U_before = (m1 * m2) / r_col_norm
   Etot_before = KE_before + U_before
   KE_after = 0.0_DP
   U_after = 0.0_DP
   ! Calculate the position and velocity of each fragment 
   do i=1, nfrag 
      v_frag_norm = sqrt(2 * f_anelastic * Etot_before / (nfrag * m_frag(i))) 
      r_frag_norm = r_col_norm * mtot / m_frag(i) 

      v_frag(:,i) =  v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:) + &
                                 (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))
      x_frag(:,i) =  r_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)  + &
                                 (sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))

      !KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
   end do
   !do i = 1, nfrag
   !   do j = i+1, nfrag
   !      delta_x(:) = x_frag(:,i) - x_frag(:,j)
   !      U_after = U_after + m_frag(i) * m_frag(j) / norm2(delta_x(:))
   !   end do
   !end do 

   !write(*,*) "SYMBA_FRAG_POS KE_before : ", KE_before
   !write(*,*) "SYMBA_FRAG_POS KE_after : ", KE_after
   !write(*,*) "SYMBA_FRAG_POS KE_ratio : ", KE_after / KE_before
   !write(*,*) "SYMBA_FRAG_POS U_before : ", U_before
   !write(*,*) "SYMBA_FRAG_POS U_after : ", U_after
   !write(*,*) "SYMBA_FRAG_POS U_ratio : ", U_after / U_before
   !write(*,*) "SYMBA_FRAG_POS E_before : ", KE_before + U_before
   !write(*,*) "SYMBA_FRAG_POS E_after : ", KE_after + U_after
   !write(*,*) "SYMBA_FRAG_POS E_ratio : ", (KE_after + U_after) / (KE_before + U_before)
!
   return 
end subroutine symba_frag_pos
