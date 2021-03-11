!**********************************************************************************************************************************
!
!  Unit Name   : symba_fragment_calculation
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  DescrIption : Compute the position of added fragments after the end of the step
!
!  Input
!    Arguments : 
!            
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : 
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL
!
!  Notes       : Adapted from Hal Levison's Swift routine util_hills.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_fragment_calculation(nmergeadd, mergeadd_list, symba_plA, plplenc_list, index_enc)

! Modules
   USE swiftest
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_helio
   USE module_symba
   USE module_swiftestalloc
   USE module_interfaces, EXCEPT_THIS_ONE => symba_fragment_calculation
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: nmergeadd, index_enc
   TYPE(symba_merger), INTENT(INOUT)                      :: mergeadd_list
   TYPE(symba_pl), INTENT(INOUT)                          :: symba_plA
   TYPE(symba_plplenc), INTENT(INOUT)                     :: plplenc_list

! Internals

   INTEGER(I4B)                                           :: frags_added, i, j, name1, name2, name_keep, name_rm
   REAL(DP)                                               :: phase_ang, r_circle, theta, v_frag_norm, v_col_norm, r_col_norm
   REAL(DP)                                               :: rhill_p1, rhill_p2, m1, m2, mass_rm, mass_keep, mtot, spin_vec_mag_frag
   REAL(DP)                                               :: r_frag, Ip_frag, r1, r2, r_keep, r_rm
   INTEGER(I4B), DIMENSION(2)                             :: idx
   REAL(DP), DIMENSION(NDIM)                              :: xh_1, xh_2, vb_1, vb_2, xh_keep, xh_rm, vb_keep, vb_rm, vh_keep, vh_rm
   REAL(DP), DIMENSION(NDIM)                              :: p_com, p_com_frag, p_com_rm, p_com_keep, p_com_1, p_com_2, p_f
   REAL(DP), DIMENSION(NDIM)                              :: v_com, v_com_frag, v_com_rm, v_com_keep, v_com_1, v_com_2, v_f
   REAL(DP), DIMENSION(NDIM)                              :: v_col_vec, v_col_unit_vec, tri_pro, tri_pro_unit_vec, rot_rm, rot_keep
   REAL(DP), DIMENSION(NDIM)                              :: vbs, delta_v, delta_p, v_cross_p, xv_rm, xv_keep, Ip_rm, Ip_keep
   REAL(DP), DIMENSION(NDIM)                              :: xv_1, xv_2, Ip_1, Ip_2, rot_1, rot_2, pv_frag, spin_hat_frag, tmp
   REAL(DP), DIMENSION(NDIM)                              :: l_orb_before, l_orb_after, l_spin_before, l_spin_after, l_spin_frag
   REAL(DP), DIMENSION(NDIM)                              :: mv_frag, mp_frag
   REAL(DP), DIMENSION(:, :), ALLOCATABLE                 :: p_frag, v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE                    :: m_frag 
   integer(I4B), save                                     :: thetashift = 0
   integer(I4B), parameter                                :: SHIFTMAX = 9
   INTEGER(I4B), DIMENSION(10)                            :: nmergeadd_frag_index

! Executable code

   ! Bring in all of the information on the two bodies involved in the collision
   idx(1) = plplenc_list%index1(index_enc)
   idx(2) = plplenc_list%index2(index_enc)
   if (symba_plA%helio%swiftest%mass(idx(2)) > symba_plA%helio%swiftest%mass(idx(1))) then
      idx(1) = plplenc_list%index2(index_enc)
      idx(2) = plplenc_list%index1(index_enc)
   end if

   vb_1 = symba_plA%helio%swiftest%vb(:, idx(1))
   vb_2 = symba_plA%helio%swiftest%vb(:, idx(2))
   xh_1 = symba_plA%helio%swiftest%xh(:, idx(1))
   xh_2 = symba_plA%helio%swiftest%xh(:, idx(2))
   m1 = symba_plA%helio%swiftest%mass(idx(1))
   m2 = symba_plA%helio%swiftest%mass(idx(2))
   rhill_p1 = symba_plA%helio%swiftest%rhill(idx(1))
   rhill_p2 = symba_plA%helio%swiftest%rhill(idx(2))
   name1 = symba_plA%helio%swiftest%name(idx(1))
   name2 = symba_plA%helio%swiftest%name(idx(2))
   r1 = symba_plA%helio%swiftest%radius(idx(1))
   r2 = symba_plA%helio%swiftest%radius(idx(2))
   Ip_1 = symba_plA%helio%swiftest%Ip(:, idx(1))
   Ip_2 = symba_plA%helio%swiftest%Ip(:, idx(2))
   rot_1 = symba_plA%helio%swiftest%rot(:, idx(1))
   rot_2 = symba_plA%helio%swiftest%rot(:, idx(2))

   nmergeadd_frag_index(:) = 0
   frags_added = 0
   vbs(:) = symba_plA%helio%swiftest%vb(:, 1)

   ! Loop through all the new bodiess in mergeadd_list and check their parents' names
   DO i = 1, nmergeadd
      ! If both of their parents' names match the two bodies we are considering in this collision
      ! then we know that this new body in mergeadd_list formed from this collision 
      IF ((mergeadd_list%name_p1(i) == name1) .AND. (mergeadd_list%name_p2(i) == name2)) THEN
         frags_added = frags_added + 1 ! Count up the fragments that formed from this collision
         nmergeadd_frag_index(frags_added) = i ! Index in mergeadd_list of the fragment from this collision
         ! If both of their parents' names match the name of the removed body in this collision
         ! then we know that this new body in mergeadd_list is the removed body from the collision 
      ELSE IF ((mergeadd_list%name_p1(i) == name1) .AND. (mergeadd_list%name_p2(i) == name1)) THEN
         ! Assign this body the position and velocity it had at the end of the step as if there was no collision
         mergeadd_list%xh(:,i) = xh_1(:)
         mergeadd_list%vh(:,i) = vb_1(:) - vbs(:)
      ELSE IF ((mergeadd_list%name_p1(i) == name2) .AND. (mergeadd_list%name_p2(i) == name2)) THEN
         mergeadd_list%xh(:,i) = xh_2(:)
         mergeadd_list%vh(:,i) = vb_2(:) - vbs(:)
      END IF
   END DO

   IF (frags_added == 0) RETURN

   ALLOCATE(m_frag(frags_added))
   ALLOCATE(p_frag(NDIM, frags_added))
   ALLOCATE(v_frag(NDIM, frags_added))

   ! Find COM of the collision
   p_com(:) = ((xh_1(:) * m1) + (xh_2(:) * m2)) / (m1 + m2)
   p_com_1(:) = xh_1 - p_com
   p_com_2(:) = xh_2 - p_com
   v_com(:) = ((vb_1(:) * m1) + (vb_2(:) * m2)) / (m1 + m2)
   v_com_1(:) = vb_1 - v_com
   v_com_2(:) = vb_2 - v_com
   delta_v(:) = vb_2(:) - vb_1(:)
   delta_p(:) = xh_2(:) - xh_1(:)

   ! Find collision velocity
   v_col_norm = NORM2(delta_v) ! collision velocity magnitude
   r_col_norm = NORM2(delta_p) ! collision velocity magnitude
   v_col_unit_vec(:) = delta_v(:) / v_col_norm ! unit vector of collision velocity (direction only)

   ! Determine the radius and angular spacing of fragments placed in a circle around the COM of the collision
   r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added))
   theta = (2 * PI) / frags_added
   ! Shifts the starting circle of fragments around so that multIple fragments generated 
   ! from a single collision in a single time step don't pile up on top of each other
   phase_ang = theta * thetashift / SHIFTMAX
   thetashift = thetashift + 1
   IF (thetashift >= shiftmax) thetashift = 0
   ! Fragment velocity is the  collision velocity with an adjustment for new distance using vis viva
   v_frag_norm = sqrt(v_col_norm**2 - 2 * (m1 + m2) * (1._DP / r_col_norm - 1._DP / r_circle))

   call util_crossproduct(p_com_1, v_com_1, xv_1)
   call util_crossproduct(p_com_2, v_com_2, xv_2)

   ! Now that all the fragment positions and velocities have been calculated, we can calculate the spins
   ! Calculate the orbital angular momentum and the spin angular momentum of the two colliding bodies before the collision
   l_orb_before = (m1 * xv_1) + (m2 * xv_2)
   l_spin_before = (Ip_1 * m1 * r1**2 * rot_1) + (Ip_2 * m2 * r2**2 * rot_2)

   m_frag(:) = 0.0_DP
   p_frag(:,:) = 0.0_DP
   v_frag(:,:) = 0.0_DP
   mv_frag(:) = 0.0_DP
   mp_frag(:) = 0.0_DP

   ! Add the masses of all fragments to m_frag and calculate the total mass of the fragments
   m_frag(1:frags_added) = mergeadd_list%mass(nmergeadd_frag_index(1):nmergeadd_frag_index(frags_added))
   mtot = sum(m_frag(1:frags_added))

   ! Calculate the trIple product
   call util_crossproduct(delta_v,delta_p,v_cross_p)
   call util_crossproduct(v_cross_p,delta_v,tri_pro)

   tri_pro_unit_vec(:) = tri_pro(:) / NORM2(tri_pro(:))
   ! Calculate the position and velocity of each fragment 
   DO i=1, frags_added ! fragment velocity (same mag for each just different direction)
      v_frag(:,i) = v_frag_norm * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)) + &
                                  ((sin(phase_ang + theta * i)) * tri_pro_unit_vec(:))
      p_frag(:,i) = r_circle    * ((cos(phase_ang + theta * i)) * v_col_unit_vec(:)) + &
                                  ((sin(phase_ang + theta * i)) * tri_pro_unit_vec)
      mv_frag(:) = mv_frag(:) + (v_frag(:,i) * m_frag(i))
      mp_frag(:) = mp_frag(:) + (p_frag(:,i) * m_frag(i))
   END DO
   ! Loop through all the fragments in this collision and add their positions and velocities to mergeadd_list
   DO i=1, frags_added
      mergeadd_list%vh(:,nmergeadd_frag_index(1)+i-1) = v_frag(:, i) - vbs(:) + v_com(:)
      mergeadd_list%xh(:,nmergeadd_frag_index(1)+i-1) = p_frag(:, i) + p_com(:)
      write(*,*) 'frag index: ',nmergeadd_frag_index(1)+i-1
   END DO

   l_orb_after(:) = 0.0_DP
   tmp(:) = 0.0_DP
   spin_vec_mag_frag = 0.0_DP
    
   ! Loop through all the fragments in this collision       
   DO i = 1, frags_added
      ! Calculate the orbital angular momentum of each fragment
      call util_crossproduct(p_frag(:,i), v_frag(:,i), pv_frag)
      ! Loop through each dimension of the orbital angular momentum 
      l_orb_after(:) = l_orb_after(:) + (m_frag(i) * pv_frag(:))
   END DO
   ! Calculate the spin angular momentum of the collisional system after collision through conservation of angular momentum
   ! AKA whatever angular momentum is lost by the orbit, is picked up by the spin
   l_spin_after = l_orb_before + l_spin_before - l_orb_after
   l_spin_frag = l_spin_after / frags_added ! The amount of spin angular momentum that each fragment will have
   spin_hat_frag = l_spin_after / (NORM2(l_spin_after)) ! The unit vector of the spin angular momentum that each fragment will have

   ! Loop through all the fragments in this collision 
   DO i = 1, frags_added
      r_frag = mergeadd_list%radius(nmergeadd_frag_index(1)+i-1)
      Ip_frag = 2.0_DP / 5.0_DP ! Because each body is a perfect sphere, the princIpal moments of inertia for each fragment will be the same
      mergeadd_list%Ip(:,nmergeadd_frag_index(1)+i-1) = Ip_frag
      ! Calculate the magnitude of the spin angular momentum of each fragment 
      spin_vec_mag_frag = norm2(l_spin_frag)  / (Ip_frag  * m_frag(i) * r_frag**2)
      ! Calculate the final spin (rotation) that each fragment will have in all three dimensions
      mergeadd_list%rot(:,nmergeadd_frag_index(1)+i-1) = spin_vec_mag_frag*spin_hat_frag(:)
      tmp(:) = tmp(:) + (spin_vec_mag_frag * spin_hat_frag(:) * Ip_frag * m_frag(i) * r_frag**2)
   END DO 

   DEALLOCATE(p_frag)
   DEALLOCATE(m_frag)
   DEALLOCATE(v_frag)
   
   RETURN

END SUBROUTINE symba_fragment_calculation