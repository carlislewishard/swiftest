subroutine symba_frag_adjust (xcom, vcom, x, v, mass, radius, L_spin, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
   !! 
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, EXCEPT_THIS_ONE => symba_frag_adjust
   implicit none

   real(DP), dimension(:),   intent(in)    :: xcom, vcom, mass, radius
   real(DP), dimension(:,:), intent(in)    :: x, v, L_spin
   real(DP), dimension(:,:), intent(in)    :: Ip_frag
   real(DP), dimension(:), intent(in)      :: m_frag, rad_frag
   real(DP), dimension(:,:), intent(out)   :: x_frag, v_frag, rot_frag

   real(DP), dimension(NDIM, 2)            :: rot, Ip
   integer(I4B)                            :: i, j, nfrag
   real(DP)                                :: mtot
   real(DP), dimension(NDIM)               :: xc, vc, x_cross_v, v_phi_unit, h_unit, v_r_unit
   real(DP), dimension(NDIM)               :: L_orb_old, L_spin_old, L_orb_new, L_spin_frag, L_residual, L_spin_new
   real(DP), dimension(NDIM)               :: mx_frag, mv_frag, COM_offset_x, COM_offset_v

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

   ! Adjust the position and velocity of the fragments as needed to align them with the original trajectory center of mass
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
end subroutine symba_frag_adjust
