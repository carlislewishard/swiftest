!**********************************************************************************************************************************
!
!  Unit Name   : symba_casehitandrun
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_casehitandrun(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, vbs, & 
   symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mres, m1, m2, rad1, rad2, xh_1, xh_2, vb_1, vb_2, mtiny)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_swiftestalloc
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casehitandrun
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc
   INTEGER(I4B), INTENT(IN)                         :: nplplenc
   INTEGER(I4B), INTENT(INOUT)                      :: plmaxname, tpmaxname, nmergeadd, nmergesub
   REAL(DP), INTENT(IN)                             :: t, dt, mtiny
   REAL(DP), INTENT(INOUT)                          :: m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres
   REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: xh_1, xh_2, vb_1, vb_2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added, j, index1_child, index2_child 
   INTEGER(I4B)                                     :: index1_parent, index2_parent, index_keep_parent, index_rm_parent
   INTEGER(I4B)                                     :: name1, name2, index_keep, index_rm, name_keep, name_rm, nstart
   REAL(DP)                                         :: mtot, msun, d_rm, m_rm, r_rm
   REAL(DP)                                         :: rhill_keep, r_circle, theta, radius1, radius2, e, q, semimajor_encounter
   REAL(DP)                                         :: m_rem, mass1, mass2, semimajor_inward
   REAL(DP)                                         :: mass_keep, mass_rm, rhill_rm
   REAL(DP)                                         :: rad_keep, rad_rm
   REAL(DP)                                         :: r_smallestcircle
   REAL(DP)                                         :: phase_ang, v_col_norm, v2esc, v2el, v2esc_circle
   REAL(DP), DIMENSION(:, :), ALLOCATABLE           :: v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE              :: m_frag
   REAL(DP), DIMENSION(NDIM)                        :: mv, xh_keep, xh_rm, vh_keep, vh_rm, xbs
   REAL(DP), DIMENSION(NDIM)                        :: vb_keep, vb_rm, tri_pro, delta_v, v_cross_p, tri_pro_unit_vec
   REAL(DP), DIMENSION(NDIM)                        :: v_com, xr, v_col_vec, v_col_unit_vec, mv_frag, v_com_frag, v_f
   INTEGER(I4B), DIMENSION(NCHILDMAX)               :: array_index1_child, array_index2_child
   INTEGER(I4B), SAVE                               :: thetashift = 0
   INTEGER(I4B), PARAMETER                          :: SHIFTMAX = 9

! Executable code

   ! Set the maximum number of fragments to be added in a Hit and Run collision (nfrag)
   nfrag = 4
   call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  

   ! Pull in the information about the two particles involved in the collision 
   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_plA%index_parent(index1)
   index2_parent = symba_plA%index_parent(index2)
   name1 = symba_plA%helio%swiftest%name(index1)
   name2 = symba_plA%helio%swiftest%name(index2)
   mass1 = symba_plA%helio%swiftest%mass(index1) ! the mass of the first particle in the collision NOT INCLUDING all it's children
   mass2 = symba_plA%helio%swiftest%mass(index2)
   radius1 = symba_plA%helio%swiftest%radius(index1)
   radius2 = symba_plA%helio%swiftest%radius(index2)
   msun = symba_plA%helio%swiftest%mass(1)
   xbs(:) = symba_plA%helio%swiftest%xb(:,1)

   ! Determine which of the two particles in the collision is larger where mass INCLUDES the mass of all their children
   IF (m2 > m1) THEN
      index_keep = index2
      index_rm = index1
      mass_keep = m2
      mass_rm = m1
      rad_keep = rad2
      rad_rm = rad1
      xh_keep = xh_2 
      xh_rm = xh_1 
      vb_keep = vb_2 
      vb_rm = vb_1 
      vh_keep = vb_rm - vbs
      vh_rm = vb_keep - vbs
      index_keep_parent = index2_parent
      index_rm_parent = index1_parent
      name_keep = name2
      name_rm = name1
   ELSE
      index_keep = index1
      index_rm = index2
      mass_keep = m1
      mass_rm = m2
      rad_keep = rad1
      rad_rm = rad2
      xh_keep = xh_1 
      xh_rm = xh_2 
      vb_keep = vb_1 
      vb_rm = vb_2 
      vh_keep = vb_keep - vbs
      vh_rm = vb_rm - vbs
      index_keep_parent = index1_parent
      index_rm_parent = index2_parent
      name_keep = name1
      name_rm = name2
   END IF

   ! Find energy pre-frag

   ! Go through the encounter list and look for particles actively encoutering in this timestep
   ! Prevent them from having further encounters in this timestep by setting status in plplenc_list to MERGED
   ! DO k = 1, nplplenc 
   !    IF ((plplenc_list%status(k) == ACTIVE) .AND. &
   !       ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
   !       (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
   !          plplenc_list%status(k) = MERGED
   !    END IF
   ! END DO

   ! Set the status of the particles in symba_plA to HIT_AND_RUN
   symba_plA%helio%swiftest%status(index1) = HIT_AND_RUN
   symba_plA%helio%swiftest%status(index2) = HIT_AND_RUN

   symba_plA%helio%swiftest%status(index1_parent) = HIT_AND_RUN
   symba_plA%helio%swiftest%status(index2_parent) = HIT_AND_RUN
   array_index1_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,index1_parent)
   array_index2_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,index2_parent)
   DO k = 1, nplplenc                                          !go through the encounter list and for particles actively encoutering, get their children
      IF (plplenc_list%status(k) == ACTIVE) THEN
         DO i = 0, symba_plA%nchild(index1_parent)
            IF (i == 0) THEN 
               index1_child = index1_parent
            ELSE
               index1_child = array_index1_child(i)
            END IF 
            DO j = 0, symba_plA%nchild(index2_parent)
               IF (j == 0) THEN
                  index2_child = index2_parent
               ELSE
                  index2_child = array_index2_child(j)
               END IF
               IF ((index1_child == plplenc_list%index1(k)) .OR. (index2_child == plplenc_list%index2(k))) THEN
                  plplenc_list%status(k) = MERGED
               ELSE IF ((index1_child == plplenc_list%index2(k)) .OR. (index2_child == plplenc_list%index1(k))) THEN
                  plplenc_list%status(k) = MERGED
               END IF
            END DO
         END DO
      END IF
   END DO

   mtot = 0.0_DP ! running total mass of new fragments
   mv(:) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd ! start of new fragments in mergeadd_list
   ! Increment around the circle for positions of fragments
   ! Calculate the positions of the new fragments in a circle of radius rhill_keep
   rhill_keep = symba_plA%helio%swiftest%rhill(index_keep_parent)
   rhill_rm = symba_plA%helio%swiftest%rhill(index_rm_parent)
   r_smallestcircle = (rhill_keep + rhill_rm) / (2.0_DP*sin(PI /2.0_DP))

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(xh_keep, vb_keep, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_smallestcircle)) THEN
      WRITE(*,*) "WARNING in symba_casehitandrun: Timestep is too large to resolve fragments."
   END IF

   ! The largest fragment = the kept parent particle
   nmergeadd = nmergeadd + 1
   mergeadd_list%status(nmergeadd) = HIT_AND_RUN
   mergeadd_list%ncomp(nmergeadd) = 2
   mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_keep)
   mergeadd_list%mass(nmergeadd) = mass_keep
   mergeadd_list%radius(nmergeadd) = rad_keep
   mergeadd_list%vh(:,nmergeadd) = vh_keep(:)
   mtot = mtot + mergeadd_list%mass(nmergeadd) 

   ! Pure Hit & Run
   IF ((mres(2) > mass_rm * 0.9_DP).OR.(mres(2) < nfrag * mtiny)) THEN
      !frags_added does NOT get incremented on in a perfect merger because then plmaxname would be plmaxname + 1
      !this screws up the naming of new fragments in subsequent disruptions or supercatastrophic disruptions or
      !imperfect hit & runs. In other words, in a hit & run, frags_added is only incremented on in imperfect 
      !hit & runs. 
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_rm)
      mergeadd_list%mass(nmergeadd) = mass_rm
      mergeadd_list%radius(nmergeadd) = rad_rm
      mergeadd_list%vh(:,nmergeadd) = vh_rm(:)
      mtot = mtot + mergeadd_list%mass(nmergeadd)

   ELSE
      m_rm = mass_rm
      r_rm = rad_rm
      d_rm = (3 * m_rm) / (4 * PI * r_rm**3)
      frags_added = frags_added + 1
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) +  1
      mergeadd_list%mass(nmergeadd) = mres(2)
      mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * d_rm))  & 
            ** (1.0_DP / 3.0_DP) 
      mtot = mtot + mergeadd_list%mass(nmergeadd)
   ! Imperfect Hit & Run
      m_rem = m1 + m2 - mass_keep - mres(2)      
      DO i = 2, nfrag 
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%status(nmergeadd) = HIT_AND_RUN
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + i
            mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
            mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * d_rm))  & 
               ** (1.0_DP / 3.0_DP) 
            mtot = mtot + mergeadd_list%mass(nmergeadd)
         END DO
   END IF

   allocate(m_frag(frags_added + 1))
   allocate(v_frag(NDIM, frags_added))

   IF (frags_added > 0) THEN

         !!!!!!!!!!!!                     DEV                      !!!!!!!!!!!!!!!! 
         r_circle = (rhill_keep + rhill_rm) / (2.0_DP*sin(PI / frags_added))
         theta = (2 * PI) / (frags_added)
         m_frag(1:frags_added + 1) = mergeadd_list%mass(nstart + 1 :nstart + frags_added + 1)
         write(*,*) "casehitandrun m_frag", m_frag
         
         mtot = sum(m_frag(1:frags_added))
         mv_frag = 0.0_DP

         ! Shifts the starting circle of fragments around so that multiple fragments generated 
         ! in from a single body in a single time step don't pile up on top of each other
         phase_ang = theta * thetashift / SHIFTMAX
         thetashift = thetashift + 1
         if (thetashift >= shiftmax) thetashift = 0

         ! Find velocity of COM
         v_com(:) = ((vb_keep(:) * mass_keep) + (vb_rm(:) * mass_rm)) / (mass_keep + mass_rm)

         ! Find Collision velocity
         xr(:) = xh_rm(:) - xh_keep(:) ! distance between particles at time of collision
         v_col_norm = NORM2(vb_rm(:) - vb_keep(:)) ! collision velocity magnitude
         v_col_vec(:) = (vb_rm(:) - vb_keep(:)) ! collision velocity vector
         v_col_unit_vec(:) = v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)

         v2esc = 2.0_DP * GC * (mass_keep+mass_rm) / (NORM2(xr(:))) ! escape velocity from COM squared
         v2esc_circle = 2.0_DP * GC * (mass_keep+mass_rm) * (1.0_DP/(NORM2(xr)) - 1.0_DP/r_circle) ! escape velocity from circle squared
         v2el = - SQRT(v2esc - v2esc_circle) ! adjusted escape velocity to account for distance from COM

         ! Calculate the triple product
         delta_v(:) = vb_rm(:) - vb_keep(:)

         call util_crossproduct(delta_v,xr,v_cross_p)
         call util_crossproduct(v_cross_p,delta_v,tri_pro)

         tri_pro_unit_vec(:) = tri_pro(:) / NORM2(tri_pro(:))

         ! Calculate the velocity magnitude and direction of each fragment 
         DO i=1, frags_added ! fragment velocity (same mag for each just different direction)
            v_frag(:,i) = ((v2el * cos(phase_ang + theta * i))*v_col_unit_vec(:)) + &
            ((v2el * sin(phase_ang + theta + i)) * tri_pro_unit_vec) + v_com(:)
            mv_frag(:) = (v_frag(:,i) * m_frag(i)) + mv_frag(:) ! rolling linear momentum of the system
         END DO

         ! Calculate the error 
         v_com_frag(:) = mv_frag(:) / mtot ! velocity of the COM of the fragments
         v_f(:) = v_com(:) - v_com_frag(:) ! velocity error between COM of collison and COM of fragments

         DO i=1, frags_added
            v_frag(:,i) = v_frag(:,i) + v_f(:) ! velocity of the fragments including the error
            mergeadd_list%vh(:, nstart + i) = v_frag(:, i) - vbs(:) ! add to mergeadd_list 
            write(*,*) "casehitandrun fragnum, vh", i, mergeadd_list%vh(:, nstart + i)
         END DO

         !!!!!!!!!!!!                     DEV                      !!!!!!!!!!!!!!!! 
   END IF

   deallocate(m_frag)
   deallocate(v_frag)

   ! Add both particles involved in the collision to mergesub_list
   call symba_merger_size_check(mergesub_list, nmergesub + 2) 
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = HIT_AND_RUN 
   mergesub_list%xh(:,nmergesub) = xh_keep
   mergesub_list%vh(:,nmergesub) = vh_keep !v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = rad1
   IF (frags_added == 0) THEN !AKA if it was a perfect merger
      mergesub_list%nadded(nmergesub) = 2
   ELSE 
      mergesub_list%nadded(nmergesub) = frags_added + 1 !the plus one is from the biggest body
   END IF
   mergesub_list%index_ps(nmergesub) = index_keep

   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = HIT_AND_RUN
   mergesub_list%xh(:,nmergesub) = xh_rm
   mergesub_list%vh(:,nmergesub) = vh_rm !v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = rad2
   IF (frags_added == 0) THEN !AKA if it was a perfect merger
      mergesub_list%nadded(nmergesub) = 2
   ELSE 
      mergesub_list%nadded(nmergesub) = frags_added + 1 !the plus one is from the biggest body
   END IF
   mergesub_list%index_ps(nmergesub) = index_rm

   WRITE(*, *) "Hit and run between particles ", name1, " and ", name2, " at time t = ",t
   IF (frags_added == 0) THEN
      WRITE(*,*) "0 fragments produced; pure hit and run."
   ELSE
      WRITE(*, *) "Particle ", name_keep, " survives; Particle ", name_rm, " is fragmented."
      WRITE(*, *) "Number of fragments added: ", (frags_added)
   END IF

   ! Update plmaxname to account for new fragments made in imperfect hit & runs
   plmaxname = max(plmaxname, tpmaxname) + frags_added
   RETURN 
END SUBROUTINE symba_casehitandrun

