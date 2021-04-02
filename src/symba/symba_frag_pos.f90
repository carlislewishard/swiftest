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
   type(symba_pl)                            :: symba_pla

   real(DP), dimension(NDIM, 2)            :: rot
   integer(I4B)                            :: i, j, nfrag
   real(DP), dimension(NDIM)               :: v_cross_x, delta_v, delta_x, L_spin_tot, xcom, vcom
   real(DP)                                :: mtot, phase_ang, theta, v_frag_norm, r_frag_norm, v_col_norm, r_col_norm
   real(DP)                                :: f_anelastic, Etot_before, KE_before, U_before, v1mag2, v2mag2, U_p1, U_p2
   real(DP)                                :: U_after, KE_spin_before, KE_spin_after, KE_after, KE_corrected, U_corrected
   real(DP), dimension(NDIM)               :: v_col_unit_vec, tri_pro, tri_pro_unit_vec, v_com, vc1, vc2
   integer(I4B), save                      :: thetashift = 0
   integer(I4B), parameter                 :: SHIFTMAX = 9
   real(DP), dimension(2)                  :: m

   nfrag = size(x_frag, 2)
   mtot = sum(mass(:))
   xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
   vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot 

   m(:) = symba_plA%helio%swiftest%mass(idx_parents(:)) ! just the parents, not including their children

   ! Now work out where to put the new bodies such that we conserve momentum and lose an appropriate amount of energy
   ! Find collision velocity
   delta_v(:) = v(:, 2) - v(:, 1)
   delta_x(:) = x(:, 2) - x(:, 1)

   v_com(:) = ((v(:, 1) * mass(1)) + (v(:, 2) * mass(2))) / mtot

   v_col_norm = norm2(delta_v(:)) ! pre-collision velocity magnitude
   r_col_norm = norm2(delta_x(:)) ! pre-collision distance 
   v_col_unit_vec(:) = delta_v(:) / v_col_norm ! unit vector of collision velocity 

   f_anelastic = 0.1_DP ! TODO: Should this be set by the user or kept as a constant?
   vc1(:) = v(:,1) - v_com(:)
   vc2(:) = v(:,2) - v_com(:)
   v1mag2 = dot_product(vc1(:), vc1(:))
   v2mag2 = dot_product(vc2(:), vc2(:))

   ! Add in kinetic energy from the orbit and from spin
   KE_before = 0.5_DP * (mass(1) * v1mag2 + mass(2) * v2mag2) 
   rot(:,1) = L_spin(:,1) / (mass(1) * Ip(3,1) * radius(1)**2)
   rot(:,2) = L_spin(:,2) / (mass(2) * Ip(3,2) * radius(2)**2)
   L_spin_tot(:) = L_spin(:,1) + L_spin(:,2)
   KE_spin_before = 0.5_DP * (dot_product(rot(:,1),L_spin(:,1)) + dot_product(rot(:,2),L_spin(:,2)))

   U_p1 = 0.0_DP
   U_p2 = 0.0_DP
   !! Go through all the bodies in the system and add up their potential energy contribution
   do i = 1, size(symba_plA%helio%swiftest%mass(:))
      !! Skip a body in symba_plA if it is one of the parent bodies
      if ((symba_plA%helio%swiftest%name(i) == symba_plA%helio%swiftest%name(idx_parents(1))) &
         .or. (symba_plA%helio%swiftest%name(i) == symba_plA%helio%swiftest%name(idx_parents(1)))) cycle
      !! Add up the potential energy contribution of each body on parent 1
      U_p1 = U_p1 - (m(1) * symba_plA%helio%swiftest%mass(i) / norm2(x(:,1) - symba_plA%helio%swiftest%xh(:,i))) 
      !! Add up the potential energy contribution of each body on parent 2
      U_p2 = U_p2 - (m(2) * symba_plA%helio%swiftest%mass(i) / norm2(x(:,2) - symba_plA%helio%swiftest%xh(:,i))) 
   end do

   U_before = - (m(1) * m(2) / r_col_norm) - U_p1 - U_p2
   Etot_before = KE_before + KE_spin_before + U_before

   ! Calculate the position of each fragment 
   ! Theta is a phase shift value that ensures that successive nearby collisions in a single step are rotated to avoid possible overlap
   theta = (2 * PI) / nfrag
   ! Shifts the starting circle of fragments around so that multiple fragments generated 
   ! from a single collision in a single time step don't pile up on top of each other
   phase_ang = theta * thetashift / SHIFTMAX
   thetashift = thetashift + 1
   IF (thetashift >= shiftmax) thetashift = 0

   ! Calculate the triple product to get the plane of the fragment distribution
   call util_crossproduct(delta_v,delta_x,v_cross_x)
   call util_crossproduct(v_cross_x,delta_v,tri_pro)

   tri_pro_unit_vec(:) = tri_pro(:) / norm2(tri_pro(:))

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
   U_after = U_before

   ! Adjust fragment positiosn so that their center of mass follows the original pre-impact trajectory
   call symba_frag_adjust (xcom, vcom, x, v, mass, radius, L_spin, Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag)
   
   ! Compute the current kinetic energy of the system so we can scale the velocity vectors to the correct value
   KE_after = 0.0_DP
   KE_spin_after = 0.0_DP
   do i = 1, nfrag
      KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
      KE_spin_after = 0.5_DP * m_frag(i) * rad_frag(i)**2 * Ip_frag(3,i) * dot_product(rot_frag(:,i), rot_frag(:,i))
   end do

   ! Adjust the fragment velocities so that they have the their total energy reduced by an amount set by the anelastic parameter


   ! Make sure we don't end up with negative energy (bound system). If so, we'll adjust the radius so that the potential energy
   ! takes up the negative part
   KE_corrected = f_anelastic * Etot_before - KE_spin_after - U_after
   if (KE_corrected < 0.0_DP) then
      U_corrected = U_after + KE_corrected
      x_frag(:,:) = x_frag(:,:) * U_corrected / U_after
      U_after = U_corrected
      KE_corrected = 0.0_DP
   end if

   v_frag(:,:) = v_frag(:,:) * sqrt(KE_corrected / KE_after)

   ! REMOVE THE FOLLOWING AFTER TESTING
   ! Calculate the new energy of the system of fragments
   KE_after = 0.0_DP
   do i = 1, nfrag
      KE_after = KE_after + 0.5_DP * m_frag(i) * dot_product(v_frag(:,i), v_frag(:,i))
   end do
   U_after = 0.0_DP
   do j = 1, nfrag - 1
      do i = j+1, nfrag
         delta_x(:) = x_frag(:,i) - x_frag(:,j)
         U_after = U_after - m_frag(i) * m_frag(j) / norm2(delta_x(:))
      end do 
   end do

   write(*,*) "SYMBA_FRAG_POS KE_before : ", KE_before + KE_spin_before
   write(*,*) "SYMBA_FRAG_POS KE_after  : ", KE_after + KE_spin_after
   write(*,*) "SYMBA_FRAG_POS KE_ratio  : ", (KE_after + KE_spin_after) / (KE_before + KE_spin_before)
   write(*,*) "SYMBA_FRAG_POS U_before  : ", U_before
   write(*,*) "SYMBA_FRAG_POS U_after   : ", U_after
   write(*,*) "SYMBA_FRAG_POS U_ratio   : ", U_after / U_before
   write(*,*) "SYMBA_FRAG_POS E_before  : ", KE_before + KE_spin_before + U_before
   write(*,*) "SYMBA_FRAG_POS E_after   : ", KE_after + KE_spin_after + U_after
   write(*,*) "SYMBA_FRAG_POS E_ratio   : ", (KE_after + KE_spin_after + U_after) / (KE_before + KE_spin_before + U_before)
   read(*,*)

   do i = 1, nfrag
      x_frag(:, i) = x_frag(:, i) + xcom(:)
      v_frag(:, i) = v_frag(:, i) + vcom(:)
   end do


   return 
end subroutine symba_frag_pos
