!**********************************************************************************************************************************
!
!  Unit Name   : symba_casedisruption
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
!  Invocation  : CALL symba_casedisruption(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, vbs, & 
   symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mres, rres, m1, m2, rad1, rad2, xh_1, xh_2, vb_1, vb_2)

! Modules
   USE swiftest
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_helio
   USE module_symba
   USE module_swiftestalloc
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casedisruption
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc
   INTEGER(I4B), INTENT(IN)                         :: nplplenc
   INTEGER(I4B), INTENT(INOUT)                      :: plmaxname, tpmaxname, nmergeadd, nmergesub
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: xh_1, xh_2, vb_1, vb_2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added, index1_child, index2_child 
   INTEGER(I4B)                                     :: index1_parent, index2_parent, j 
   INTEGER(I4B)                                     :: name1, name2, nstart
   REAL(DP)                                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   REAL(DP)                                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2, r_smallestcircle
   REAL(DP)                                         :: m_rem, mass1, mass2
   REAL(DP)                                         :: phase_ang, v_col_norm, v2esc, v2el, v2esc_circle
   REAL(DP), DIMENSION(:, :), ALLOCATABLE           :: v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE              :: m_frag
   REAL(DP), DIMENSION(NDIM)                        :: mv, xbs, vh_1, vh_2
   REAL(DP), DIMENSION(NDIM)                        :: v_com, xr, v_col_vec, v_col_unit_vec, mv_frag, v_com_frag, v_f
   INTEGER(I4B), DIMENSION(NCHILDMAX)               :: array_index1_child, array_index2_child
   INTEGER(I4B), SAVE                               :: thetashift = 0
   INTEGER(I4B), PARAMETER                          :: SHIFTMAX = 9

! Executable code

   ! Set the maximum number of fragments to be added in a Disruption collision (nfrag)
   nfrag = 5 
   call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  

   ! Pull in the information about the two particles involved in the collision 
   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_plA%index_parent(index1)
   index2_parent = symba_plA%index_parent(index2)
   name1 = symba_plA%helio%swiftest%name(index1)
   name2 = symba_plA%helio%swiftest%name(index2)
   mass1 = symba_plA%helio%swiftest%mass(index1) ! the mass of the first particle in the collision NOT INCLUDING all its children
   mass2 = symba_plA%helio%swiftest%mass(index2)
   radius1 = symba_plA%helio%swiftest%radius(index1)
   radius2 = symba_plA%helio%swiftest%radius(index2)
   msun = symba_plA%helio%swiftest%mass(1)
   xbs(:) = symba_plA%helio%swiftest%xb(:,1)
   vh_1(:) = vb_1(:) - vbs(:)
   vh_2(:) = vb_2(:) - vbs(:)

   WRITE(*, *) "Disruption between particles ", name1, " and ", name2, " at time t = ",t

   ! Set the status of the particles in symba_plA to DISRUPTION
   symba_plA%helio%swiftest%status(index1) = DISRUPTION
   symba_plA%helio%swiftest%status(index2) = DISRUPTION
   symba_plA%helio%swiftest%status(index1_parent) = DISRUPTION
   symba_plA%helio%swiftest%status(index2_parent) = DISRUPTION
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
   
   rhill_p1 = symba_plA%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_plA%helio%swiftest%rhill(index2_parent)
   r_smallestcircle = (rhill_p1 + rhill_p2) / (2.0_DP*sin(PI /2.0_DP))

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(xh_1, vb_1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_smallestcircle)) THEN
      WRITE(*,*) "WARNING in symba_casedisruption: Timestep is too large to resolve fragments."
   END IF
   ! If not, continue through all possible fragments to be added
   mtot = 0.0_DP ! running total mass of new fragments
   mv(:) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd ! start of new fragments in mergeadd_list

   d_p1 = (3 * m1) / (4 * PI * rad1**3)
   d_p2 = (3 * m2) / (4 * PI * rad2**3)
   avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)

   frags_added = frags_added + 1
   nmergeadd = nmergeadd + 1
   mergeadd_list%status(nmergeadd) = DISRUPTION
   mergeadd_list%ncomp(nmergeadd) = 2
   mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + frags_added 
   mergeadd_list%mass(nmergeadd) = mres(1)
   mergeadd_list%radius(nmergeadd) = rres(1)
   mtot = mtot + mergeadd_list%mass(nmergeadd)

   IF ((mres(2) > (1.0_DP / 3.0_DP)*mres(1))) THEN !DM to JP and CW: What is the purpose of this line?
      ! frags_added is the actual number of fragments added to the simulation vs nfrag which is the total possible
      frags_added = frags_added + 1
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = DISRUPTION
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + frags_added
      mergeadd_list%mass(nmergeadd) = mres(2)
      mergeadd_list%radius(nmergeadd) = rres(2)
      mtot = mtot + mergeadd_list%mass(nmergeadd)
      DO i = 3, nfrag
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%status(nmergeadd) = DISRUPTION
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + frags_added
         m_rem = (m1 + m2) - (mres(1) + mres(2))
         mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 2) 
         mtot = mtot + mergeadd_list%mass(nmergeadd) 
         if (i == nfrag) then
            ! If there is any residual mass left at the end, put it in the last body
            m_rem = (m1 + m2) - mtot
            mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
         end if
         mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * avg_d))  & 
            ** (1.0_DP / 3.0_DP)
      END DO                           

   ELSE   
      DO i = 2, nfrag
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%status(nmergeadd) = DISRUPTION
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + frags_added
        
         m_rem = (m1 + m2) - mres(1)
         mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
         mtot = mtot + mergeadd_list%mass(nmergeadd)
         if (i == nfrag) then
         ! If there is any residual mass left at the end, put it in the last body
            m_rem = (m1 + m2) - mtot
            mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
         end if
         mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * avg_d))  & 
            ** (1.0_DP / 3.0_DP)  
      END DO
   END IF

   !!!!!!!!!!!!                     DEV                      !!!!!!!!!!!!!!!! 

   ! Calculate the positions of the new fragments in a circle with a radius large enough to space
   ! all fragments apart by a distance of rhill_p1 + rhill_p2
   r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added)) !((2.0_DP * rhill_p1 + 2.0_DP * rhill_p2) / (2.0_DP * sin(PI / frags_added))) 
   theta = (2 * PI) / frags_added

   ALLOCATE(m_frag(frags_added))
   ALLOCATE(v_frag(NDIM, frags_added))
   m_frag(1:frags_added) = mergeadd_list%mass(nstart + 1 :nstart + frags_added)

   mtot = sum(m_frag(1:frags_added))
   m_rem = (m1 + m2) - mtot
   mv_frag = 0.0_DP

   ! Shifts the starting circle of fragments around so that multiple fragments generated 
   ! in from a single body in a single time step don't pile up on top of each other
   phase_ang = theta * thetashift / SHIFTMAX
   thetashift = thetashift + 1
   if (thetashift >= shiftmax) thetashift = 0

   ! Find velocity of COM
   v_com(:) = ((vb_1(:) * m1) + (vb_2(:) * m2)) / (m1 + m2)

   ! Find Collision velocity
   xr(:) = xh_2(:) - xh_1(:) ! distance between particles at time of collision
   v_col_norm = NORM2(vb_2(:) - vb_1(:)) ! collision velocity magnitude
   v_col_vec(:) = (vb_2(:) - vb_1(:)) ! collision velocity vector
   v_col_unit_vec(:) = v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)

   v2esc = 2.0_DP * GC * (m1+m2) / (NORM2(xr(:))) ! escape velocity from COM squared
   v2esc_circle = 2.0_DP * GC * (m1+m2) * (1.0_DP/(NORM2(xr)) - 1.0_DP/r_circle) ! escape velocity from circle squared
   v2el = - SQRT(v2esc - v2esc_circle) ! adjusted escape velocity to account for distance from COM

   ! Calculate the velocity magnitude and direction of each fragment 
   DO i=1, frags_added
      v_frag(:,i) = ((v2el * cos(phase_ang + theta * i))*v_col_unit_vec(:)) + v_com(:) ! fragment velocity (same mag for each just different direction)
      mv_frag(:) = (v_frag(:,i) * m_frag(i)) + mv_frag(:) ! rolling linear momentum of the system
   END DO

   ! Calculate the error 
   v_com_frag(:) = mv_frag(:) / mtot ! velocity of the COM of the fragments
   v_f(:) = v_com(:) - v_com_frag(:) ! velocity error between COM of collison and COM of fragments

   ! Calculate the final velocity of the fragments and add them to mergeadd_list
   DO i=1, frags_added
      v_frag(:,i) = v_frag(:,i) + v_f(:) ! velocity of the fragments including the error
      mergeadd_list%vh(:, nstart + i) = v_frag(:, i) - vbs(:) ! add to mergeadd_list 
   END DO 

   deallocate(m_frag)
   deallocate(v_frag)

   !!!!!!!!!!!!                     DEV                      !!!!!!!!!!!!!!!!   

   ! Add both particles involved in the collision to mergesub_list
   call symba_merger_size_check(mergesub_list, nmergesub + 2)  
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = DISRUPTION
   mergesub_list%xh(:,nmergesub) = xh_1(:)!x1(:)
   mergesub_list%vh(:,nmergesub) = vh_1(:)!v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   mergesub_list%nadded(nmergesub) = frags_added
   mergesub_list%index_ps(nmergesub) = index1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = DISRUPTION
   mergesub_list%xh(:,nmergesub) = xh_2(:)!x2(:)
   mergesub_list%vh(:,nmergesub) = vh_2(:)!v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = radius2
   mergesub_list%nadded(nmergesub) = frags_added
   mergesub_list%index_ps(nmergesub) = index2

   WRITE(*, *) "Number of fragments added: ", frags_added

   plmaxname = max(plmaxname, tpmaxname) + frags_added
   RETURN 
END SUBROUTINE symba_casedisruption