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
SUBROUTINE symba_frag_pos(nmergeadd_step, nmergesub_step, nmergeadd, nmergesub, mergeadd_list, mergesub_list)

! Modules
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_interfaces, EXCEPT_THIS_ONE => symba_mom
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: nmergeadd_step, nmergesub_step, nmergeadd, nmergesub
   TYPE(symba_merger), INTENT(INOUT)                      :: mergeadd_list, mergesub_list

! Internals

   INTEGER(I4B)                                           :: count_enc, count_frag, numenc, nmergeadd_start, nmergesub_start
   REAL(DP)                                               :: phase_ang, r_circle, rhill_p1, rhill_p2, m1, m2, r1, r2, v_com_norm
   REAL(DP)                                               :: m_frag_tot
   REAL(DP), DIMENSION(NDIM)                              :: p_com, v_col_vec, v_col_unit_vec, mp_frag, p_com_frag, p_f
   REAL(DP), DIMENSION(NDIM)                              :: xh_1, xh_2, vh_1, vh_2
   REAL(DP), DIMENSION(:, :), ALLOCATABLE                 :: p_frag
   REAL(DP), DIMENSION(:), ALLOCATABLE                    :: m_frag
   integer(I4B), save                                     :: thetashift = 0
   integer(I4B), parameter                                :: SHIFTMAX = 9


! Executable code

   numenc = nmergesub_step / 2 !number of encounters this step
   nmergesub_start = nmergesub - nmergesub_step + 1 !where the particles subtracted in this step are located in mergesub_list
   nmergeadd_start = nmergeadd - nmergeadd_step + 1 !where the particles added in this step are located in mergeadd_list

   count_enc = 0 !counter for the number of encountering bodies in this timestep used to increment on mergesub_list

   DO i = 1, numenc
      ! First particle in encounter pair
      xh_1(:) = mergesub_list%xh(nmergesub_start + count_enc)
      vh_1(:) = mergesub_list%vh(nmergesub_start + count_enc)
      m1 = mergesub_list%mass(nmergesub_start + count_enc)
      r1 = mergesub_list%radius(nmergesub_start + count_enc)
      rhill_p1 = 
      ! Second particle in encounter pair
      xh_2(:) = mergesub_list%xh(nmergesub_start + count_enc + 1)
      vh_2(:) = mergesub_list%vh(nmergesub_start + count_enc + 1)
      m2 = mergesub_list%mass(nmergesub_start + count_enc + 1)
      r2 = mergesub_list%radius(nmergesub_start + count_enc + 1)
      rhill_p2 = 

      frags_added = mergesub_list%nadded(nmergesub_start + count_enc)

      ALLOCATE(m_frag(frags_added))
      ALLOCATE(p_frag(NDIM, frags_added))
   
      ! Calculate the positions of the new fragments in a circle with a radius large enough to space
      ! all fragments apart by a distance of rhill_p1 + rhill_p2
      r_circle = (rhill_p1 + rhill_p2) / (2 * sin(PI / frags_added)) !((2.0_DP * rhill_p1 + 2.0_DP * rhill_p2) / (2.0_DP * sin(PI / frags_added))) 
      theta = (2 * PI) / frags_added

      ! Shifts the starting circle of fragments around so that multiple fragments generated in from a single body in a single time step 
      ! don't pile up on top of each other
      phase_ang = theta * thetashift / SHIFTMAX
      thetashift = thetashift + 1
      if (thetashift >= shiftmax) thetashift = 0

      ! Find COM
      p_com(:) = ((xh_1(:) * m1) + (xh_2(:) * m2)) / (m1 + m2)

      ! Find Collision velocity
      v_col_norm = NORM2(vb_2(:) - vb_1(:)) ! collision velocity magnitude
      v_col_vec(:) = (vb_2(:) - vb_1(:)) ! collision velocity vector
      v_col_unit_vec(:) = v_col_vec(:) / v_col_norm ! unit vector of collision velocity (direction only)

      mp_frag = 0.0_DP

      count_frag = 0 !counter for the number of fragments added in this timestep used to increment on mergeadd_list

      DO i=1, frags_added
         m_frag(i) = mergeadd_list%mass(nmergeadd_start + count_frag + i - 1)
         p_frag(:,i) = (- r_circle  * cos(phase_ang + theta * i))*v_col_unit_vec(:) + p_com(:)
         mp_frag = (p_frag(:,i) * m_frag(i)) + mp_frag(:)
      END DO

      m_frag_tot = SUM(m_frag(:))
      p_com_frag(:) = mp_frag(:) / m_frag_tot
      p_f(:) =  p_com(:) - p_com_frag(:)

      DO i=1, frags_added
         p_frag(:,i) = p_frag(:,i) + p_f(:)
      END DO 

      count_frag = count_frag + frags_added

      DEALLOCATE(p_frag)
      DEALLOCATE(m_frag)

      count_enc = count_enc + 2
   END DO 

   RETURN

END SUBROUTINE symba_frag_pos