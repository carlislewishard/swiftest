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

   real(DP), dimension(:),   intent(in)      :: xcom, vcom, mass, radius
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin
   real(DP), dimension(:,:), intent(in)      :: Ip_frag
   real(DP), dimension(:), intent(in)        :: m_frag, rad_frag
   real(DP), dimension(:,:), intent(out)     :: x_frag, v_frag, rot_frag

   real(DP), dimension(NDIM, 2)            :: rot, Ip
   integer(I4B)                            :: i, j, nfrag
   real(DP)                                :: mtot
   real(DP), dimension(NDIM)               :: xc, vc, x_cross_v, vc1, vc2, delta_x
   real(DP), dimension(NDIM)               :: L_orb_old, L_spin_old, L_spin_new, L_orb_new, L_spin_frag, L_spin_tot
   real(DP), dimension(NDIM)               :: mx_frag, mv_frag, COM_offset_x, COM_offset_v
   real(DP)                                :: Etot_before, KE_before, U_before, v1mag2, v2mag2
   real(DP)                                :: U_after, KE_spin_after, KE_after

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

   ! Calculate the spin angular momentum of the collisional system after collision through conservation of angular momentum
   ! AKA whatever angular momentum is lost by the orbit, is picked up by the spin of all the fragments
   L_orb_new(:) = 0.0_DP
   do i = 1, nfrag
      call utiL_crossproduct(x_frag(:, i), v_frag(:, i), x_cross_v(:))
      L_orb_new(:) = L_orb_new(:) + m_frag(i) * x_cross_v(:)
   end do

   L_spin_new(:) = L_orb_old(:) + L_spin_old(:) - L_orb_new(:)
   ! Now divide up the angular momentum equally between the various bodies by mass
   do i = 1, nfrag
      L_spin_frag(:) = L_spin_new(:) * m_frag(i) / mtot
      rot_frag(:,i) = L_spin_frag(:) / (Ip_frag(3, i) * m_frag(i) * rad_frag(i)**2)
   end do

   return 
end subroutine symba_frag_adjust
