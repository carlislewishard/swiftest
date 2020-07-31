!**********************************************************************************************************************************
!
!  Unit Name   : symba_casesupercatastrophic
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
!  Invocation  : CALL symba_casesupercatastrophic(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
   symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mres, rres, m1, m2, rad1, rad2, xh_1, xh_2, vb_1, vb_2)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_swiftestalloc
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casesupercatastrophic
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc
   INTEGER(I4B), INTENT(IN)                         :: nplplenc
   INTEGER(I4B), INTENT(INOUT)                      :: plmaxname, tpmaxname, nmergeadd, nmergesub
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: xh_1, xh_2, vb_1, vb_2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added, j
   INTEGER(I4B)                                     :: index1_parent, index2_parent, index1_child, index2_child 
   INTEGER(I4B)                                     :: name1, name2, nstart
   REAL(DP)                                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   REAL(DP)                                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2, r_smallestcircle
   REAL(DP)                                         :: mass1, mass2, enew, eold
   REAL(DP)                                         :: m1m2_10
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, xbs, vh_1, vh_2
   REAL(DP), DIMENSION(:, :), ALLOCATABLE           :: x_frag, v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE              :: m_frag
   INTEGER(I4B), DIMENSION(NCHILDMAX)   :: array_index1_child, array_index2_child



! Executable code
   
   ! Set the maximum number of fragments to be added in a Supercatastrophic Disruption collision (nfrag)
   nfrag = 10
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
   vh_1(:) = vb_1(:) - vbs(:) !symba_plA%helio%swiftest%vh(:,index1)
   vh_2(:) = vb_2(:) - vbs(:) !symba_plA%helio%swiftest%vh(:,index2)

   ! Find energy pre-frag
   eold = 0.5_DP * (m1 * DOT_PRODUCT(vb_1(:), vb_1(:)) + m2 * DOT_PRODUCT(vb_2(:), vb_2(:)))
   xr(:) = xh_2(:) - xh_1(:)
   eold = eold - m1 * m2 / NORM2(xr(:)) 

   WRITE(*, *) "Supercatastrophic disruption between particles ", name1, " and ", name2, " at time t = ",t

   ! Go through the encounter list and look for particles actively encoutering in this timestep
   ! Prevent them from having further encounters in this timestep by setting status in plplenc_list to MERGED
   ! DO k = 1, nplplenc 
   !    IF ((plplenc_list%status(k) == ACTIVE) .AND. &
   !       ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
   !       (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
   !          plplenc_list%status(k) = MERGED
   !    END IF
   ! END DO

   ! Set the status of the particles in symba_plA to DISRUPTION
   symba_plA%helio%swiftest%status(index1) = SUPERCATASTROPHIC
   symba_plA%helio%swiftest%status(index2) = SUPERCATASTROPHIC
   symba_plA%helio%swiftest%status(index1_parent) = SUPERCATASTROPHIC
   symba_plA%helio%swiftest%status(index2_parent) = SUPERCATASTROPHIC
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

   ! Calculate the positions of the new fragments in a circle with a radius large enough to space
   ! all fragments apart by a distance of rhill_p1 + rhill_p2
   rhill_p1 = symba_plA%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_plA%helio%swiftest%rhill(index2_parent)
   r_smallestcircle = (rhill_p1 + rhill_p2) / (2.0_DP*sin(PI /2.0_DP))

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(xh_1, vb_1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_smallestcircle)) THEN
      WRITE(*,*) "WARNING in symba_casesupercatastrophic: Timestep is too large to resolve fragments."
   END IF
   ! If not, continue through all possible fragments to be added
      mtot = 0.0_DP ! running total mass of new fragments
      mv(:) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
      frags_added = 0 ! running total number of new fragments
      m1m2_10 = 0.1_DP * (m1 + m2) ! one tenth the total initial mass of the system used to check the size of the fragments
      nstart = nmergeadd

      d_p1 = (3 * m1) / (4 * PI * rad1**3)
      d_p2 = (3 * m2) / (4 * PI * rad2**3)
      avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)

         ! If we are adding the first and largest fragment (lr), check to see if its mass is SMALLER than one tenth the total
         ! mass of the system aka if it is too small to resolve. If so, add a fragment with a mass of one tenth the total mass 
         ! of the system and calculate its radius.
      IF ((mres(1) < m1m2_10)) THEN
         DO i = 1, nfrag
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + i
            mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%mass(nmergeadd) = m1m2_10
            mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * avg_d))  & 
                              ** (1.0_DP / 3.0_DP)
            mtot = mtot + mergeadd_list%mass(nmergeadd) 
         END DO 
      ELSE
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + 1
         mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%mass(nmergeadd) = mres(1)
         mergeadd_list%radius(nmergeadd) = rres(1)
         mtot = mtot + mergeadd_list%mass(nmergeadd) 
         ! Fragments creation 
         DO i = 2, nfrag
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%name(nmergeadd) = max(plmaxname, tpmaxname) + i
            mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%mass(nmergeadd) = (m1 + m2 - mres(1)) / (nfrag - 1)
            mergeadd_list%radius(nmergeadd) = ((3 * mergeadd_list%mass(nmergeadd)) / (4 * PI * avg_d))  & 
                  ** (1.0_DP / 3.0_DP)  
            mtot = mtot + mergeadd_list%mass(nmergeadd)
         END DO
      END IF 

   r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added))
   theta = (2 * PI) / frags_added

   ALLOCATE(m_frag(frags_added))

   m_frag(1:frags_added) = mergeadd_list%mass(nstart + 1:nstart + frags_added)

   ALLOCATE(x_frag(NDIM, frags_added))
   ALLOCATE(v_frag(NDIM, frags_added))
   CALL symba_mom(m1, xh_1, vb_1, m2, xh_2, vb_2, frags_added, m_frag, r_circle, theta, x_frag, v_frag)

   DO i=1, frags_added

      mergeadd_list%xh(:,nstart + i) = x_frag(:, i)! - xbs(1)!x_frag
      mergeadd_list%vh(:,nstart + i) = v_frag(:, i) - vbs(:)!vx_frag

         ! Tracking linear momentum. 
      mv(:) = mv(:) + (mergeadd_list%mass(nstart + i) * mergeadd_list%vh(:,nstart + i))
   END DO
   deallocate(m_frag)
   deallocate(x_frag)
   deallocate(v_frag)

   ! Add both particles involved in the collision to mergesub_list
   call symba_merger_size_check(mergesub_list, nmergesub + 2)  
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = SUPERCATASTROPHIC
   mergesub_list%xh(:,nmergesub) = xh_1(:)!x1(:)
   mergesub_list%vh(:,nmergesub) = vh_1(:)!v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   mergesub_list%nadded(nmergesub) = frags_added
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = SUPERCATASTROPHIC
   mergesub_list%xh(:,nmergesub) = xh_2(:)!x2(:)
   mergesub_list%vh(:,nmergesub) = vh_2(:)!v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = radius2
   mergesub_list%nadded(nmergesub) = frags_added

   WRITE(*, *) "Number of fragments added: ", frags_added
   ! Calculate energy after frag                                                                           
   vnew(:) = mv(:) / mtot    ! COM of new fragments                               
   enew = 0.5_DP * mtot * DOT_PRODUCT(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew
   ! Update plmaxname to account for new fragments
   plmaxname = max(plmaxname, tpmaxname) + frags_added

   RETURN 
END SUBROUTINE symba_casesupercatastrophic

