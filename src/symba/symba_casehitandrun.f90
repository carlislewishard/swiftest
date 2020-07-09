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
SUBROUTINE symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
   symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casehitandrun
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
   INTEGER(I4B), INTENT(IN)                         :: nplplenc
   INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, fragmax
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: x1, x2, v1, v2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added
   INTEGER(I4B)                                     :: index1_parent, index2_parent, index_keep_parent, index_rm_parent
   INTEGER(I4B)                                     :: name1, name2, index_keep, index_rm, name_keep, name_rm, nstart
   real(DP)                                         :: first_add_vz, second_add_vz, first_add_pz, second_add_pz
   real(DP)                                         :: first_add_name, second_add_name
   REAL(DP)                                         :: mtot, msun, d_rm, m_rm, r_rm, x_rm, y_rm, z_rm, vx_rm, vy_rm, vz_rm 
   REAL(DP)                                         :: rhill_keep, r_circle, theta, radius1, radius2, e, q, semimajor_encounter
   REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold, semimajor_inward, A, B, v_col
   REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com, mass_keep, mass_rm, rhill_rm
   REAL(DP)                                         :: rad_keep, rad_rm
   REAL(DP)                                         :: r_smallestcircle
   REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE     :: x_frag, v_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE        :: m_frag
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, xh_keep, xh_rm, vh_keep, vh_rm, xbs

! Executable code

   ! Set the maximum number of fragments to be added in a Hit and Run collision (nfrag)
   nfrag = 4
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
      xh_keep = x2
      xh_rm = x1
      vh_keep = v2
      vh_rm = v1
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
      xh_keep = x1
      xh_rm = x2
      vh_keep = v1
      vh_rm = v2
      index_keep_parent = index1_parent
      index_rm_parent = index2_parent
      name_keep = name1
      name_rm = name2
   END IF

   ! Find COM
   x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! Find Collision velocity
   v_col = NORM2(v2(:) - v1(:))

   ! Find energy pre-frag
   eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
   xr(:) = x2(:) - x1(:)
   eold = eold - (m1*m2/(SQRT(DOT_PRODUCT(xr(:), xr(:)))))

   ! Go through the encounter list and look for particles actively encoutering in this timestep
   ! Prevent them from having further encounters in this timestep by setting status in plplenc_list to MERGED
   DO k = 1, nplplenc 
      IF ((plplenc_list%status(k) == ACTIVE) .AND. &
         ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
         (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
            plplenc_list%status(k) = MERGED
      END IF
   END DO

   ! Set the status of the particles in symba_plA to HIT_AND_RUN
   symba_plA%helio%swiftest%status(index1) = HIT_AND_RUN
   symba_plA%helio%swiftest%status(index2) = HIT_AND_RUN


   mtot = 0.0_DP ! running total mass of new fragments
   mv(1) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   mv(2) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   mv(3) = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd + 1 ! start of new fragments in mergeadd_list
   ! Increment around the circle for positions of fragments
   ! Calculate the positions of the new fragments in a circle of radius rhill_keep
   rhill_keep = symba_plA%helio%swiftest%rhill(index_keep_parent)
   rhill_rm = symba_plA%helio%swiftest%rhill(index_rm_parent)
   r_smallestcircle = (rhill_keep + rhill_rm) / (2.0_DP*sin(PI /2.0_DP))

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
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
   mergeadd_list%xh(:,nmergeadd) = xh_keep(:)
   mergeadd_list%vh(:,nmergeadd) = vh_keep(:)
   mtot = mtot + mergeadd_list%mass(nmergeadd) 


   ! Pure Hit & Run
   IF (mres(2) > mass_rm * 0.9_DP) THEN
      !frags_added does NOT get incremented on in a perfect merger because then fragmax would be fragmax + 1
      !this screws up the naming of new fragments in subsequent disruptions or supercatastrophic disruptions or
      !imperfect hit & runs. In other words, in a hit & run, frags_added is only incremented on in imperfect 
      !hit & runs. 
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_rm)
      mergeadd_list%mass(nmergeadd) = mass_rm
      mergeadd_list%radius(nmergeadd) = rad_rm
      mergeadd_list%xh(:,nmergeadd) = xh_rm(:)
      mergeadd_list%vh(:,nmergeadd) = vh_rm(:)
      mtot = mtot + mergeadd_list%mass(nmergeadd)
   ELSE
      m_rm = mass_rm
      r_rm = rad_rm
      d_rm = (3.0_DP * m_rm) / (4.0_DP * PI * (r_rm ** 3.0_DP))
      frags_added = frags_added + 1
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + 1
      mergeadd_list%mass(nmergeadd) = mres(2)
      mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * d_rm))  & 
            ** (1.0_DP / 3.0_DP) 
      mtot = mtot + mergeadd_list%mass(nmergeadd)
   ! Imperfect Hit & Run       
      DO i = 2, nfrag 
         m_rem = m_rm - mres(2)
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%status(nmergeadd) = HIT_AND_RUN
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
         mergeadd_list%mass(nmergeadd) = m_rem / (nfrag) 
         mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * d_rm))  & 
            ** (1.0_DP / 3.0_DP) 
         mtot = mtot + mergeadd_list%mass(nmergeadd)
      END DO
   END IF

   IF (frags_added > 0) THEN
         r_circle = (rhill_keep + rhill_rm) / (2.0_DP*sin(PI / frags_added))
         theta = (2.0_DP * PI) / (frags_added)
         ALLOCATE(m_frag(frags_added))
         m_frag(1:frags_added) = mergeadd_list%mass(nstart + 1 :nstart + 1 + frags_added)

         ALLOCATE(x_frag(NDIM, frags_added))
         ALLOCATE(v_frag(NDIM, frags_added))
         CALL util_mom(0.0_DP, xh_keep+xbs, vh_keep, mass_rm, xh_rm+xbs, vh_rm, & 
            frags_added, nstart, m_frag, r_circle, theta, x_frag, v_frag)
         DO i=1, frags_added

            mergeadd_list%xh(1,nstart + i) = x_frag(1, i) -xbs(1) !x_frag
            mergeadd_list%xh(2,nstart + i) = x_frag(2, i) -xbs(2) !y_frag
            mergeadd_list%xh(3,nstart + i) = x_frag(3, i) -xbs(3) !z_frag                                                   
            mergeadd_list%vh(1,nstart + i) = v_frag(1, i) -vbs(1) !vx_frag
            mergeadd_list%vh(2,nstart + i) = v_frag(2, i) -vbs(2) !vy_frag
            mergeadd_list%vh(3,nstart + i) = v_frag(3, i) -vbs(1)  !vz_frag

         ! Tracking linear momentum. 
            mv(:) = mv(:) + (mergeadd_list%mass(nstart + i) * mergeadd_list%vh(:,nstart + i))
         END DO
         deallocate(m_frag)
         deallocate(x_frag)
         deallocate(v_frag)
   END IF

   ! Add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = HIT_AND_RUN 
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = rad1
   IF (frags_added == 0) THEN !AKA if it was a perfect merger
      !You must have nadded for a pure hit & run be equal to 1 so it does the discard correctly
      mergesub_list%nadded(nmergesub) = 1
   ELSE 
      mergesub_list%nadded(nmergesub) = frags_added
   END IF
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = HIT_AND_RUN
   mergesub_list%xh(:,nmergesub) = x2(:)
   mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = rad2
   IF (frags_added == 0) THEN !AKA if it was a perfect merger
      mergesub_list%nadded(nmergesub) = 2
   ELSE 
      mergesub_list%nadded(nmergesub) = frags_added
   END IF

   WRITE(*, *) "Hit and run between particles ", name1, " and ", name2, " at time t = ",t
   IF (frags_added == 0) THEN
      WRITE(*,*) "0 fragments produced; pure hit and run."
   ELSE
      WRITE(*, *) "Particle ", name_keep, " survives; Particle ", name_rm, " is fragmented."
      WRITE(*, *) "Number of fragments added: ", (frags_added)
   END IF

   first_add_name = mergeadd_list%name(1)
   second_add_name = mergeadd_list%name(2)

   first_add_pz = mergeadd_list%xh(3,1)
   second_add_pz = mergeadd_list%xh(3,2)

   first_add_vz = mergeadd_list%vh(3,1)
   second_add_vz = mergeadd_list%vh(3,2)
   
   ! Calculate energy after frag                                                                           
   vnew(:) = mv(:) / mtot    ! COM of new fragments                               
   enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew
   ! Update fragmax to account for new fragments made in imperfect hit & runs
   fragmax = fragmax + frags_added
   RETURN 
END SUBROUTINE symba_casehitandrun

