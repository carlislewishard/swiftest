!**********************************************************************************************************************************
!
!  Unit Name   : symba_frag_pos
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
SUBROUTINE symba_frag_pos(nmergeadd_step, nmergesub_step, nmergeadd, nmergesub, mergeadd_list, mergesub_list, symba_plA, npl)

! Modules
   USE swiftest
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_helio
   USE module_symba
   USE module_swiftestalloc
   USE module_interfaces, EXCEPT_THIS_ONE => symba_frag_pos
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: nmergeadd_step, nmergesub_step, nmergeadd, nmergesub, npl
   TYPE(symba_merger), INTENT(INOUT)                      :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                          :: symba_plA

! Internals

   INTEGER(I4B)                                           :: count_enc, count_frag, numenc, nmergeadd_start, nmergesub_start, i, j, k
   INTEGER(I4B)                                           :: frags_added
   REAL(DP)                                               :: phase_ang, r_circle, rhill_p1, rhill_p2, m1, m2, r1, r2, v_col_norm
   REAL(DP)                                               :: m_frag_tot, theta,spin_vec_mag_frag
   REAL(DP), DIMENSION(NDIM)                              :: p_com, v_col_vec, v_col_unit_vec, mp_frag, p_com_frag, p_f, tri_pro
   REAL(DP), DIMENSION(NDIM)                              :: xvh_1, xvh_2, pv_frag, spin_hat_frag
   REAL(DP), DIMENSION(NDIM)                              :: xh_1, xh_2, vh_1, vh_2, vbs, vb_1, vb_2, delta_v, delta_p, v_cross_p
   REAL(DP), DIMENSION(NDIM)                              :: tri_pro_unit_vec, vh_1_end, vh_2_end, xh_rm, IP_1, IP_2, rot_1, rot_2
   REAL(DP), DIMENSION(NDIM)                              :: l_orb_before, l_orb_after, l_spin_before, l_spin_after, l_spin_frag
   REAL(DP), DIMENSION(:, :), ALLOCATABLE                 :: p_frag, v_frag, IP_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE                    :: m_frag
   integer(I4B), save                                     :: thetashift = 0
   integer(I4B), parameter                                :: SHIFTMAX = 9
   REAL(DP)                                               :: ke, pe, te, msys, Ltot_now, Ltot_after, te_orig, te_error, name_1, name_2, mass, radius, IP
   REAL(DP), DIMENSION(NDIM)                              :: htot, x, v,  rot, h, dx, htot_now, htot_after, v_f, v_com_frag, v_com, mv_frag 

! Executable code
   write(*,*) "enter symba_frag_pos"
   call symba_energy(npl, symba_plA%helio%swiftest, 0.0_DP, 0.0_DP, ke, pe, te_orig, htot, msys)
   Ltot_now = norm2(htot)
   htot_now = htot 
   write(*,*) "htot_now", htot 
   write(*,*) "Ltot_now", Ltot_now
   write(*,*) "te_orig = ", te_orig 
   numenc = nmergesub_step / 2 !number of encounters this step
   nmergesub_start = nmergesub - nmergesub_step + 1 !where the particles subtracted in this step are located in mergesub_list
   nmergeadd_start = nmergeadd - nmergeadd_step + 1 !where the particles added in this step are located in mergeadd_list

   count_enc = 0 !counter for the number of encountering bodies in this timestep used to increment on mergesub_list

   xh_1(:) = 0.0_DP
   vh_1_end(:) = 0.0_DP
   rhill_p1 = 0.0_DP

   xh_2(:) = 0.0_DP
   vh_2_end(:) = 0.0_DP
   rhill_p2 = 0.0_DP

   

   DO i = 1, numenc
      ! First particle in encounter pair
      DO j = 1, npl !loop through all the planets in symba_plA
         ! If the name of the planet in symba_plA matches the name of the planet in mergesub_list
         ! then use the position of the planet in symba_plA aka at the end of the step
         IF (symba_plA%helio%swiftest%name(j) == mergesub_list%name(nmergesub_start + count_enc)) THEN
            xh_1(:) = symba_plA%helio%swiftest%xh(:,j)
            vh_1_end(:) = symba_plA%helio%swiftest%vh(:,j)
            name_1 = symba_plA%helio%swiftest%name(j)
         END IF
      END DO
      vh_1(:) = mergesub_list%vh(:,nmergesub_start + count_enc)
      m1 = mergesub_list%mass(nmergesub_start + count_enc)
      r1 = mergesub_list%radius(nmergesub_start + count_enc) 
      rhill_p1 = symba_plA%helio%swiftest%rhill(nmergesub_start + count_enc)
      IP_1(:) = symba_plA%helio%swiftest%Ip(:,nmergesub_start + count_enc)
      rot_1(:) = symba_plA%helio%swiftest%rot(:,nmergesub_start + count_enc)
      write(*,*) "rot_1 = ", rot_1 
      write(*,*) "IP_1 = ", IP_1 
      ! Second particle in encounter pair
      DO j = 1, npl !loop through all the planets in symba_plA
         ! If the name of the planet in symba_plA matches the name of the planet in mergesub_list
         ! then use the position of the planet in symba_plA aka at the end of the step
         IF (symba_plA%helio%swiftest%name(j) == mergesub_list%name(nmergesub_start + count_enc + 1)) THEN
            xh_2(:) = symba_plA%helio%swiftest%xh(:,j)
            vh_2_end(:) = symba_plA%helio%swiftest%vh(:,j)
            name_2 = symba_plA%helio%swiftest%name(j)
         END IF
      END DO
      vh_2(:) = mergesub_list%vh(:,nmergesub_start + count_enc + 1)
      m2 = mergesub_list%mass(nmergesub_start + count_enc + 1)
      r2 = mergesub_list%radius(nmergesub_start + count_enc + 1)
      rhill_p2 = symba_plA%helio%swiftest%rhill(nmergesub_start + count_enc + 1)
      IP_2(:) = symba_plA%helio%swiftest%Ip(:,nmergesub_start + count_enc + 1)
      rot_2(:) = symba_plA%helio%swiftest%rot(:,nmergesub_start + count_enc + 1)
      write(*,*) "rot_2 = ", rot_2 
      write(*,*) "IP_2 = ", IP_2 
      frags_added = mergesub_list%nadded(nmergesub_start + count_enc)
 
      IF (frags_added > 1) THEN !if this is not a perfect merger

         ALLOCATE(m_frag(frags_added))
         ALLOCATE(p_frag(NDIM, frags_added))
         ALLOCATE(v_frag(NDIM, frags_added))
         ALLOCATE(IP_frag(NDIM, frags_added))
         m_frag(:) = 0.0_DP
         p_frag(:,:) = 0.0_DP
         IP_frag(:,:) = 0.0_DP
         mv_frag(:) = 0.0_DP

         !Calculate the positions of the new fragments in a circle with a radius large enough to space
         ! all fragments apart by a distance of rhill_p1 + rhill_p2
         IF ((mergeadd_list%status(nmergeadd_start) == HIT_AND_RUN) .and. (frags_added > 2)) THEN !this is an imperfect hit and run
            r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / (frags_added - 1))) 
            theta = (2 * PI) / (frags_added - 1)
         ELSE !this is everything else
            r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added)) 
            theta = (2 * PI) / frags_added
         END IF 

         ! Shifts the starting circle of fragments around so that multiple fragments generated in from a single body in a single time step 
         ! don't pile up on top of each other
         phase_ang = theta * thetashift / SHIFTMAX
         thetashift = thetashift + 1
         if (thetashift >= shiftmax) thetashift = 0

         ! Find COM
         p_com(:) = ((xh_1(:) * m1) + (xh_2(:) * m2)) / (m1 + m2)

         ! Find Collision velocity
         v_col_norm = NORM2(vh_2(:) - vh_1(:)) ! collision velocity magnitude
         v_col_vec(:) = (vh_2(:) - vh_1(:)) ! collision velocity vector
         v_col_unit_vec(:) = v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)

         ! Calculate the triple product
         vbs(:) = symba_plA%helio%swiftest%vb(:, 1)

         vb_1(:) = vh_1(:) + vbs(:)
         vb_2(:) = vh_2(:) + vbs(:)

         delta_v(:) = vb_2(:) - vb_1(:)
         delta_p(:) = xh_2(:) - xh_1(:)

         call util_crossproduct(delta_v,delta_p,v_cross_p)
         call util_crossproduct(v_cross_p,delta_v,tri_pro)

         tri_pro_unit_vec(:) = tri_pro(:) / NORM2(tri_pro(:))

         mp_frag = 0.0_DP

         count_frag = 0 !counter for the number of fragments added in this timestep used to increment on mergeadd_list
         
         IF ((mergeadd_list%status(nmergeadd_start) == HIT_AND_RUN) .and. (frags_added > 2)) THEN !this is an imperfect hit and run
            DEALLOCATE(m_frag)
            ALLOCATE(m_frag(frags_added + 1))
            DEALLOCATE(v_frag)
            ALLOCATE(v_frag(NDIM, frags_added + 1))

            IF (m2 > m1) THEN
               xh_rm = xh_1 
            ELSE
               xh_rm = xh_2 
            END IF

            p_com(:) = xh_rm

            DO j=1, frags_added
               m_frag(j) = mergeadd_list%mass(nmergeadd_start + count_frag + j - 1)
               v_frag(:,j) = mergeadd_list%vh(:,nmergeadd_start + count_frag + j - 1)
               p_frag(:,j) = ((- r_circle  * cos(phase_ang + theta * j)) * v_col_unit_vec(:)) + &
               ((- r_circle * sin(phase_ang + theta * j)) * tri_pro_unit_vec) + p_com(:)
               mp_frag = (p_frag(:,j) * m_frag(j)) + mp_frag(:)
            END DO
         ELSE
            DO j=1, frags_added
               m_frag(j) = mergeadd_list%mass(nmergeadd_start + count_frag + j - 1)
               v_frag(:,j) = mergeadd_list%vh(:,nmergeadd_start + count_frag + j - 1)
               p_frag(:,j) = ((- r_circle  * cos(phase_ang + theta * j)) * v_col_unit_vec(:)) + &
               ((- r_circle * sin(phase_ang + theta * j)) * tri_pro_unit_vec) + p_com(:)
               mp_frag = (p_frag(:,j) * m_frag(j)) + mp_frag(:)
               mv_frag = mv_frag + (v_frag(:,j) +vbs )*m_frag(j)  !!!! UNSURE 
            END DO
         END IF 

         m_frag_tot = SUM(m_frag(:))
         p_com_frag(:) = mp_frag(:) / m_frag_tot
         p_f(:) =  p_com(:) - p_com_frag(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNSURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Calculate the error 
         v_com(:) = ((vb_1(:) * m1) + (vb_2(:) * m2)) / (m1 + m2)
         v_com_frag(:) = mv_frag(:) / m_frag_tot ! velocity of the COM of the fragments
         v_f(:) = v_com(:) - v_com_frag(:) ! velocity error between COM of collison and COM of fragments
         ! Calculate the final velocity of the fragments and add them to mergeadd_list
         DO j=1, frags_added
            v_frag(:,j) = v_frag(:,j) + v_f(:) ! velocity of the fragments including the error
            mergeadd_list%vh(:, nmergeadd_start + count_frag + j - 1) = v_frag(:, j) ! add to mergeadd_list 
         END DO 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UNSURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ((mergeadd_list%status(nmergeadd_start) == HIT_AND_RUN) .and. (frags_added > 2)) THEN !this is an imperfect hit and run
            DO j=2, frags_added
               p_frag(:,j) = p_frag(:,j) + p_f(:)
               mergeadd_list%xh(:, nmergeadd_start + count_frag + j - 1) = p_frag(:, j)
            END DO 
            mergeadd_list%xh(:, nmergeadd_start + count_frag) = xh_1(:)
            mergeadd_list%vh(:, nmergeadd_start + count_frag) = vh_1_end(:)
         ELSE IF ((mergeadd_list%status(nmergeadd_start) == HIT_AND_RUN) .and. (frags_added == 2)) THEN !this is a perfect hit and run
            mergeadd_list%xh(:, nmergeadd_start + count_frag) = xh_1(:)
            mergeadd_list%xh(:, nmergeadd_start + count_frag + 1) = xh_2(:)
            mergeadd_list%vh(:, nmergeadd_start + count_frag) = vh_1_end(:)
            mergeadd_list%vh(:, nmergeadd_start + count_frag + 1) = vh_2_end(:)
         ELSE
            DO j=1, frags_added
               p_frag(:,j) = p_frag(:,j) + p_f(:)
               mergeadd_list%xh(:, nmergeadd_start + count_frag + j ) = p_frag(:, j)
            END DO 
         END IF

         !########################################################## DEV ###############################################################
         call util_crossproduct(xh_1, vh_1, xvh_1)
         call util_crossproduct(xh_2, vh_2, xvh_2)
         l_orb_before = (m1 * xvh_1) + (m2 * xvh_2)
         l_spin_before = (IP_1 * m1 * r1**2 * rot_1) + (IP_2 * m2 * r2**2 * rot_2)
         
         l_orb_after(:) = 0.0_DP
         DO j = 1, frags_added
            call util_crossproduct(p_frag(:,j), v_frag(:,j), pv_frag)
            DO k = 1, NDIM
               l_orb_after(k) = l_orb_after(k) + (m_frag(j) * pv_frag(k))
            END DO
         END DO

         l_spin_after = l_orb_before + l_spin_before - l_orb_after
         l_spin_frag = l_spin_after / frags_added
         spin_vec_mag_frag = 0.0_DP
         spin_hat_frag = l_spin_after / (NORM2(l_spin_after))
         DO j = 1, frags_added
            IP_frag(:,j) = (2.0_DP / 5.0_DP) 
            mergeadd_list%IP(:, nmergeadd_start + count_frag + j - 1) = IP_frag(:,j)
            spin_vec_mag_frag = norm2(l_spin_frag)  / (IP_frag(3,j)  * m_frag(j) * mergeadd_list%radius(nmergeadd_start + count_frag + j - 1)**2)
            mergeadd_list%rot(:, nmergeadd_start + count_frag + j - 1) = spin_vec_mag_frag*spin_hat_frag(:)
         END DO 
         !########################################################## DEV ################################################################
         !##########################################################CALC ###############################################################
         write(*,*) "start first loop"
         DO k = 1, npl 
            if ((symba_plA%helio%swiftest%name(k) /= name_2) .OR. (symba_plA%helio%swiftest%name(k) /= name_1)) then 
               write(*,*) "k particle = ", k 
               x(:) = symba_plA%helio%swiftest%xb(:, k)
               v(:) = symba_plA%helio%swiftest%vb(:, k)
               h(1) =  (x(2) * v(3) - x(3) * v(2))
               h(2) =  (x(3) * v(1) - x(1) * v(3))
               h(3) =  (x(1) * v(2) - x(2) * v(1))
               IP = symba_plA%helio%swiftest%Ip(3,k)
               rot = symba_plA%helio%swiftest%rot(:,k)
               radius = symba_plA%helio%swiftest%radius(k)
               htot(:) = symba_plA%helio%swiftest%mass(k) * h(:) + htot(:) + symba_plA%helio%swiftest%mass(k)*IP*rot*radius**2
               ke = ke + 0.5_DP * symba_plA%helio%swiftest%mass(k) * dot_product(v, v)
            end if 
         END DO 
         write(*,*) "first loop done"
         DO j = 1, frags_added
            x(:) = mergeadd_list%xh(:, nmergeadd_start + count_frag + j - 1)
            v(:) = mergeadd_list%vh(:, nmergeadd_start + count_frag + j - 1)
            h(1) =  (x(2) * v(3) - x(3) * v(2))
            h(2) =  (x(3) * v(1) - x(1) * v(3))
            h(3) =  (x(1) * v(2) - x(2) * v(1))
            mass = mergeadd_list%mass(nmergeadd_start + count_frag + j - 1)
            IP  = mergeadd_list%Ip(3, nmergeadd_start + count_frag + j - 1)
            rot = mergeadd_list%rot(:, nmergeadd_start + count_frag + j - 1)
            radius = mergeadd_list%radius( nmergeadd_start + count_frag + j - 1)
            htot(:) = htot(:) + mass * h(:) + mass*IP*rot*radius**2
         END DO 
         write(*,*) "second loop done"
         pe = 0.0_DP
         do k = 1, npl - 1
            if ((symba_plA%helio%swiftest%name(k) /= name_2) .OR. (symba_plA%helio%swiftest%name(k) /= name_1))then 
               do j = k + 1, npl
                  if ((symba_plA%helio%swiftest%name(j) /= name_2) .OR. (symba_plA%helio%swiftest%name(j) /= name_1)) then 
                     dx(:) = symba_plA%helio%swiftest%xh(:, j) - symba_plA%helio%swiftest%xh(:, k)
                     radius = norm2(dx(:)) 
                     pe = pe - symba_plA%helio%swiftest%mass(k) * symba_plA%helio%swiftest%mass(j) / radius
                  end if 
               end do 
               do j = 1, frags_added
                  dx(:) = mergeadd_list%xh(:, nmergeadd_start + count_frag + j - 1) - symba_plA%helio%swiftest%xb(:, k) 
                  pe = pe - symba_plA%helio%swiftest%mass(k) * mergeadd_list%mass( nmergeadd_start + count_frag + j - 1) / norm2(dx(:)) 
               end do
            end if 
         end do
         write(*,*) "third loop done"
         do j = 1, frags_added-1
            do k = j+1, frags_added
               dx(:) = mergeadd_list%xh(:, nmergeadd_start + count_frag + k - 1) - mergeadd_list%xh(:, nmergeadd_start + count_frag + j - 1)
               pe = pe - mergeadd_list%mass( nmergeadd_start + count_frag + k - 1) * mergeadd_list%mass( nmergeadd_start + count_frag + j - 1) / norm2(dx(:)) 
            END DO 
         END DO 
         write(*,*) "final loop done"
         te = ke + pe
         Ltot_after = norm2(htot)
         htot_after = htot 
         write(*,*) "htot_after", htot 
         write(*,*) "htot_after - htot_now", htot_after - htot_now
         write(*,*) "Ltot_after", (Ltot_after - Ltot_now) / Ltot_now
         write(*,*) "te = ", te 
         !##########################################################CALC ###############################################################
         count_frag = count_frag + frags_added

         DEALLOCATE(p_frag)
         DEALLOCATE(m_frag)
         DEALLOCATE(v_frag)
         DEALLOCATE(IP_frag)

      END IF

      count_enc = count_enc + 2

   END DO 

   RETURN

END SUBROUTINE symba_frag_pos