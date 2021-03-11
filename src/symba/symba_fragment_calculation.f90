!**********************************************************************************************************************************
!
!  Unit Name   : symba_fragment_calculation
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Compute the position of added fragments after the end of the step
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
   REAL(DP)                                               :: phase_ang, r_circle, theta, v2esc_circle, v2esc, v2el, v_col_norm 
   REAL(DP)                                               :: rhill_p1, rhill_p2, m1, m2, mass_rm, mass_keep, mtot, spin_vec_mag_frag
   REAL(DP)                                               :: r_frag, ip_frag, r1, r2, r_keep, r_rm
   INTEGER(I4B), DIMENSION(2)                             :: idx
   REAL(DP), DIMENSION(NDIM)                              :: xh_1, xh_2, vb_1, vb_2, xh_keep, xh_rm, vb_keep, vb_rm, vh_keep, vh_rm
   REAL(DP), DIMENSION(NDIM)                              :: p_com, p_com_frag, p_com_rm, p_com_keep, p_com_1, p_com_2, p_f
   REAL(DP), DIMENSION(NDIM)                              :: v_com, v_com_frag, v_com_rm, v_com_keep, v_com_1, v_com_2, v_f
   REAL(DP), DIMENSION(NDIM)                              :: v_col_vec, v_col_unit_vec, tri_pro, tri_pro_unit_vec, rot_rm, rot_keep
   REAL(DP), DIMENSION(NDIM)                              :: vbs, delta_v, delta_p, v_cross_p, xv_rm, xv_keep, ip_rm, ip_keep
   REAL(DP), DIMENSION(NDIM)                              :: xv_1, xv_2, ip_1, ip_2, rot_1, rot_2, pv_frag, spin_hat_frag, tmp
   REAL(DP), DIMENSION(NDIM)                              :: l_orb_before, l_orb_after, l_spin_before, l_spin_after, l_spin_frag
   REAL(DP), DIMENSION(:, :), ALLOCATABLE                 :: p_frag, v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE                    :: m_frag, mv_frag, mp_frag
   integer(I4B), save                                     :: thetashift = 0
   integer(I4B), parameter                                :: SHIFTMAX = 9
   INTEGER(I4B), DIMENSION(10)                            :: nmergeadd_frag_index

! Executable code

   ! Bring in all of the information on the two bodies involved in the collision
   idx(1) = plplenc_list%index1(index_enc)
   idx(2) = plplenc_list%index2(index_enc)

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
   ip_1 = symba_plA%helio%swiftest%ip(:, idx(1))
   ip_2 = symba_plA%helio%swiftest%ip(:, idx(2))
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
      END IF
   END DO

   ALLOCATE(m_frag(frags_added))
   ALLOCATE(p_frag(NDIM, frags_added))
   ALLOCATE(v_frag(NDIM, frags_added))
   ALLOCATE(mp_frag(NDIM))
   ALLOCATE(mv_frag(NDIM))

   m_frag(:) = 0.0_DP
   p_frag(:,:) = 0.0_DP
   v_frag(:,:) = 0.0_DP
   mv_frag(:) = 0.0_DP
   mp_frag(:) = 0.0_DP

   ! If we formed fragments in this collision AKA if it is NOT a pure hit and run
   IF (frags_added > 0) THEN
      ! Determine the radius and angular spacing of fragments placed in a circle around the COM of the collision
      r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added))
      theta = (2 * PI) / frags_added
      ! Shifts the starting circle of fragments around so that multiple fragments generated 
      ! from a single collision in a single time step don't pile up on top of each other
      phase_ang = theta * thetashift / SHIFTMAX
      thetashift = thetashift + 1
      IF (thetashift >= shiftmax) thetashift = 0
      ! Add the masses of all fragments to m_frag and calculate the total mass of the fragments
      m_frag(1:frags_added) = mergeadd_list%mass(nmergeadd_frag_index(1):nmergeadd_frag_index(frags_added))
      mtot = sum(m_frag(1:frags_added))
   END IF

   ! Calculate the COM of the fragments and the collision velocity if they are formed in a supercatastrophic or disruptive collision
   IF ((mergeadd_list%status(nmergeadd_frag_index(1)) == SUPERCATASTROPHIC) .or. (mergeadd_list%status(nmergeadd_frag_index(1)) == DISRUPTION)) THEN

      ! Find COM
      p_com(:) = ((xh_1(:) * m1) + (xh_2(:) * m2)) / (m1 + m2)
      p_com_1(:) = xh_1 - p_com
      p_com_2(:) = xh_2 - p_com
      v_com(:) = ((vb_1(:) * m1) + (vb_2(:) * m2)) / (m1 + m2)
      v_com_1(:) = vb_1 - v_com
      v_com_2(:) = vb_2 - v_com
      delta_v(:) = v_com_2(:) - v_com_1(:)
      delta_p(:) = p_com_2(:) - p_com_1(:)

      ! Find collision velocity
      v_col_norm = NORM2(v_com_2(:) - v_com_1(:)) ! collision velocity magnitude
      v_col_vec(:) = (v_com_2(:) - v_com_1(:)) ! collision velocity vector
      v_col_unit_vec(:) = v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)

      v2esc = 2.0_DP * GC * (m1+m2) / (NORM2(delta_p(:))) ! escape velocity from COM squared
      v2esc_circle = 2.0_DP * GC * (m1+m2) * (1.0_DP/(NORM2(delta_p)) - 1.0_DP/r_circle) ! escape velocity from circle squared
      v2el = - SQRT(v2esc - v2esc_circle) ! adjusted escape velocity to account for distance from COM

      call util_crossproduct(p_com_1, v_com_1, xv_1)
      call util_crossproduct(p_com_2, v_com_2, xv_2)

      ! Now that all the fragment positions and velocities have been calculated, we can calculate the spins
      ! Calculate the orbital angular momentum and the spin angular momentum of the two colliding bodies before the collision
      l_orb_before = (m1 * xv_1) + (m2 * xv_2)
      l_spin_before = (ip_1 * m1 * r1**2 * rot_1) + (ip_2 * m2 * r2**2 * rot_2)

   ! Calculate the COM of the fragments and the collision velocity if they are formed in a hit and run collision
   ELSE IF (mergeadd_list%status(nmergeadd_frag_index(1)) == HIT_AND_RUN) THEN
      ! Determine with mass is kept (the larger) and which is fragmented (the smaller)
      IF (m2 > m1) THEN
         mass_rm = m1
         mass_keep = m2
         xh_keep = xh_2 
         xh_rm = xh_1 
         vb_keep = vb_2 
         vb_rm = vb_1 
         vh_keep = vb_rm - vbs
         vh_rm = vb_keep - vbs
         name_keep = name2
         name_rm = name1
         ip_keep = ip_2
         ip_rm = ip_1
         rot_keep = rot_2
         rot_rm = rot_1
         r_keep = r2
         r_rm = r1
      ELSE
         mass_rm = m2
         mass_keep = m1
         xh_keep = xh_1 
         xh_rm = xh_2 
         vb_keep = vb_1 
         vb_rm = vb_2 
         vh_keep = vb_keep - vbs
         vh_rm = vb_rm - vbs
         name_keep = name1
         name_rm = name2
         ip_keep = ip_1
         ip_rm = ip_2
         rot_keep = rot_1
         rot_rm = rot_2
         r_keep = r1
         r_rm = r2
      END IF
      ! Loop through the planets in mergeadd_list and check their parents' names
      DO i = 1, nmergeadd
         ! If both of their parents' names match the name of the kept body in this collision
         ! then we know that this new body in mergeadd_list is the kept body from the collision 
         IF ((mergeadd_list%name_p1(i) == name_keep) .AND. (mergeadd_list%name_p2(i) == name_keep)) THEN
            ! Assign this body the position and velocity it had at the end of the step as if there was no collision
            mergeadd_list%xh(:,i) = xh_keep(:)
            mergeadd_list%vh(:,i) = vh_keep(:)
         END IF
      END DO
      ! If there are fragments formed in this hit and run collision AKA if it is an imperfect hit and run
      IF (frags_added > 0) THEN
         ! Find COM
         p_com(:) = xh_rm
         p_com_rm(:) = xh_rm - p_com
         p_com_keep(:) = xh_keep - p_com

         v_com(:) = vb_rm
         v_com_rm(:) = vb_rm - v_com
         v_com_keep(:) = vb_keep - v_com

         call util_crossproduct(p_com_rm, v_com_rm, xv_rm)
         call util_crossproduct(p_com_keep, v_com_keep, xv_keep)

         delta_v(:) = vb_keep(:) - vb_rm(:)
         delta_p(:) = xh_keep(:) - xh_rm(:)
   
         ! Find Collision velocity
         v_col_norm = NORM2(v_com_keep(:) - v_com_rm(:)) ! collision velocity magnitude
         v_col_unit_vec(:) =  delta_v(:) / v_col_norm !v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)
         v2esc = 2.0_DP * GC * mass_rm / (NORM2(delta_p(:))) ! escape velocity from COM squared
         v2esc_circle = 2.0_DP * GC * mass_rm * (1.0_DP/(NORM2(delta_p)) - 1.0_DP/r_circle) ! escape velocity from circle squared
         v2el = - SQRT(v2esc - v2esc_circle) ! adjusted escape velocity to account for distance from COM

         ! Now that all the fragment positions and velocities have been calculated, we can calculate the spins
         ! Calculate the orbital angular momentum and the spin angular momentum of the two colliding bodies before the collision
         l_orb_before = (mass_rm * xv_rm) + (mass_keep * xv_keep)
         l_spin_before = (ip_rm * mass_rm * r_rm**2 * rot_rm) + (ip_keep * mass_keep * r_keep**2 * rot_keep)

      ! If we did not form fragments in this collision AKA if it is a pure hit and run
      ELSE IF (frags_added == 0) THEN
         ! Loop through the planets in mergeadd_list and check their parents' names
         DO i = 1, nmergeadd
            ! If both of their parents' names match the name of the removed body in this collision
            ! then we know that this new body in mergeadd_list is the removed body from the collision 
            IF ((mergeadd_list%name_p1(i) == name_rm) .AND. (mergeadd_list%name_p2(i) == name_rm)) THEN
               ! Assign this body the position and velocity it had at the end of the step as if there was no collision
               mergeadd_list%xh(:,i) = xh_rm(:)
               mergeadd_list%vh(:,i) = vh_rm(:)
            END IF
         END DO
      ELSE
         ! If there is some situation in which the code gets into this loop with <0 fragments
         write(*,*) "ERROR!!! In symba_fragment_calculation, hit and run with <0 fragments!!!"
      END IF
   ELSE
      ! If there is some situation in which the code gets into this loop with no case
      write(*,*) "ERROR!!! In symba_fragment_calculation, no case selected!!!"
   END IF
   ! If we formed fragments in this collision AKA if it is NOT a pure hit and run
   IF (frags_added > 0) THEN
      ! Calculate the triple product
      call util_crossproduct(delta_v,delta_p,v_cross_p)
      call util_crossproduct(v_cross_p,delta_v,tri_pro)

      tri_pro_unit_vec(:) = tri_pro(:) / NORM2(tri_pro(:))
      ! Calculate the position and velocity of each fragment 
      DO i=1, frags_added ! fragment velocity (same mag for each just different direction)
         v_frag(:,i) = ((v2el * cos(phase_ang + theta * i))*v_col_unit_vec(:)) + &
            ((v2el * sin(phase_ang + theta + i)) * tri_pro_unit_vec)
         p_frag(:,i) = ((- r_circle  * cos(phase_ang + theta * i)) * v_col_unit_vec(:)) + &
            ((- r_circle * sin(phase_ang + theta * i)) * tri_pro_unit_vec)
         mv_frag(:) = mv_frag(:) + (v_frag(:,i) * m_frag(i))
         mp_frag(:) = mp_frag(:) + (p_frag(:,i) * m_frag(i))
      END DO
      ! Loop through all the fragments in this collision and add their positions and velocities to mergeadd_list
      DO i=1, frags_added
         mergeadd_list%vh(:,nmergeadd_frag_index(1)+i-1) = v_frag(:, i) - vbs(:) + v_com(:)
         mergeadd_list%xh(:,nmergeadd_frag_index(1)+i-1) = p_frag(:, i) + p_com(:)
      END DO
   END IF 

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
      ip_frag = 2.0_DP / 5.0_DP ! Because each body is a perfect sphere, the principal moments of inertia for each fragment will be the same
      mergeadd_list%ip(:,nmergeadd_frag_index(1)+i-1) = ip_frag
      ! Calculate the magnitude of the spin angular momentum of each fragment 
      spin_vec_mag_frag = norm2(l_spin_frag)  / (ip_frag  * m_frag(i) * r_frag**2)
      ! Calculate the final spin (rotation) that each fragment will have in all three dimensions
      mergeadd_list%rot(:,nmergeadd_frag_index(1)+i-1) = spin_vec_mag_frag*spin_hat_frag(:)
      tmp(:) = tmp(:) + (spin_vec_mag_frag * spin_hat_frag(:) * ip_frag * m_frag(i) * r_frag**2)
   END DO 

   DEALLOCATE(p_frag)
   DEALLOCATE(m_frag)
   DEALLOCATE(v_frag)
   DEALLOCATE(mp_frag)
   DEALLOCATE(mv_frag)
   
   RETURN

END SUBROUTINE symba_fragment_calculation