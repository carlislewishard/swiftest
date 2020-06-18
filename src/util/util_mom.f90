!**********************************************************************************************************************************
!
!  Unit Name   : util_mom
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Compute Hill sphere radii of planets
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_(npl, swifter_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_hills.f
!
!**********************************************************************************************************************************
SUBROUTINE util_mom(m1, xb1, vb1, m2, xb2, vb2, frags_added, nstart, m_frag, r_circle, theta, p_frag, vel_frag)

!x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag)

! Modules
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_interfaces, EXCEPT_THIS_ONE => util_mom
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: frags_added, nstart
   REAL(DP), INTENT(IN)                                   :: m1, m2, r_circle, theta
   REAL(DP), DIMENSION(:), INTENT(IN)                     :: xb1, vb1, xb2, vb2
   REAL(DP), DIMENSION(:), INTENT(IN)                     :: m_frag
   REAL(DP), DIMENSION(:, :), INTENT(INOUT)               :: p_frag, vel_frag

! Internals

   INTEGER(I4B)                                           :: i 
   REAL(DP)                                               :: x_com, y_com, z_com, vx_com, vy_com, vz_com, v_col, A
   REAL(DP)                                               :: linmom_before, angmom_after, linmom_after, DL
   REAL(DP), DIMENSION(NDIM)                              :: veclinmom_after, xbvb1, xbvb2, xv_frag, vecangmom_after
   REAL(DP)                                               :: mx_frag, my_frag, mz_frag, mvx_frag, mvy_frag, mvz_frag
   REAL(DP)                                               :: x_com_frag, y_com_frag, z_com_frag, vx_com_frag, vy_com_frag
   REAL(DP)                                               :: vz_com_frag, p_frag_check, v_frag_check, B, m_frag_tot
   REAL(DP), DIMENSION(NDIM)                              :: xr, l, kk, p, v_f, x_f, angmom_frag, angmom_fragi, angmom_com_frag
   REAL(DP), DIMENSION(NDIM)                              :: angmom_f, angmom_before

! Executable code

     linmom_before = NORM2(m1*vb1(:) + m2*vb2(:))
     call util_crossproduct(xb1,vb1,xbvb1)
     call util_crossproduct(xb2,vb2,xbvb2)
     angmom_before = (m1*xbvb1+m2*xbvb2)
     WRITE(*,*) "angmom_before =", NORM2(angmom_before)
        ! Find COM
   x_com = ((xb1(1) * m1) + (xb2(1) * m2)) / (m1 + m2)
   y_com = ((xb1(2) * m1) + (xb2(2) * m2)) / (m1 + m2)
   z_com = ((xb1(3) * m1) + (xb2(3) * m2)) / (m1 + m2)

   vx_com = ((vb1(1) * m1) + (vb2(1) * m2)) / (m1 + m2)
   vy_com = ((vb1(2) * m1) + (vb2(2) * m2)) / (m1 + m2)
   vz_com = ((vb1(3) * m1) + (vb2(3) * m2)) / (m1 + m2)

   ! Find Collision velocity
      v_col = NORM2(vb2(:) - vb1(:))

     xr(:) = xb2(:) - xb1(:)
     l(:) = (vb2(:) - vb1(:)) 
     call util_crossproduct(l(:), xr(:), kk(:))
     call util_crossproduct(kk(:), l(:), p(:))
     kk(:) = kk(:) / NORM2(kk(:))
     l(:) = l(:) / NORM2(l(:))
     p(:) = p(:) / NORM2(p(:)) 

     veclinmom_after(1) = 0.0_DP
     veclinmom_after(2) = 0.0_DP
     veclinmom_after(3) = 0.0_DP
     vecangmom_after(1) = 0.0_DP
     vecangmom_after(2) = 0.0_DP
     vecangmom_after(3) = 0.0_DP

     mx_frag = 0.0_DP
     my_frag = 0.0_DP
     mz_frag = 0.0_DP
     mvx_frag = 0.0_DP
     mvy_frag = 0.0_DP
     mvz_frag = 0.0_DP
     p_frag_check = 0.0_DP
     v_frag_check = 0.0_DP
     angmom_frag(1) = 0.0_DP
     angmom_frag(2) = 0.0_DP
     angmom_frag(3) = 0.0_DP
     A = - v_col
     B = r_circle
     DO i=1, frags_added
          p_frag(1,i) = (- B  * cos(theta * i))*l(1) + (- B  * sin(theta * i))*p(1) + x_com
          p_frag(2,i) = (- B  * cos(theta * i))*l(2) + (- B  * sin(theta * i))*p(2) + y_com
          p_frag(3,i) = (- B  * cos(theta * i))*l(3) + (- B  * sin(theta * i))*p(3) + z_com
          WRITE(*,*) "**** fragment number = ", i, " ****"
          WRITE(*,*) "p_fragx = ", p_frag(1,i)
          WRITE(*,*) "p_fragy = ", p_frag(2,i)
          WRITE(*,*) "p_fragz = ", p_frag(3,i)

          p_frag_check = - (B * cos(theta * i)) + p_frag_check
          WRITE(*,*) "pfragcheck(i)", - (r_circle * cos(theta * i)) 
          vel_frag(1,i) = (((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1)))  + vx_com
          vel_frag(2,i) = (((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2)))  + vy_com
          vel_frag(3,i) = (((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3)))  + vz_com
          WRITE(*,*) "vel_fragx = ", vel_frag(1,i)
          WRITE(*,*) "vel_fragy = ", vel_frag(2,i)
          WRITE(*,*) "vel_fragz = ", vel_frag(3,i)
          v_frag_check = (A * cos(theta * i)) + v_frag_check
          WRITE(*,*) "vfragcheck(i)", A* cos(theta * i)
          mx_frag = (p_frag(1,i) * m_frag(i)) + mx_frag
          my_frag = (p_frag(2,i) * m_frag(i)) + my_frag
          mz_frag = (p_frag(3,i) * m_frag(i)) + mz_frag

          mvx_frag = (vel_frag(1,i) * m_frag(i)) + mvx_frag
          mvy_frag = (vel_frag(2,i) * m_frag(i)) + mvy_frag
          mvz_frag = (vel_frag(3,i) * m_frag(i)) + mvz_frag

          call util_crossproduct(p_frag(:,i), vel_frag(:,i), angmom_fragi(:)) 
          angmom_frag(:) = angmom_frag(:) + angmom_fragi(:)
     END DO
     m_frag_tot = SUM(m_frag(:))
     x_com_frag = mx_frag / m_frag_tot
     y_com_frag = my_frag / m_frag_tot
     z_com_frag = mz_frag / m_frag_tot

     vx_com_frag = mvx_frag / m_frag_tot
     vy_com_frag = mvy_frag / m_frag_tot
     vz_com_frag = mvz_frag / m_frag_tot

     x_f(1) =  x_com - x_com_frag
     x_f(2) =  y_com - y_com_frag
     x_f(3) =  z_com - z_com_frag

     v_f(1) = vx_com - vx_com_frag
     v_f(2) = vy_com - vy_com_frag
     v_f(3) = vz_com - vz_com_frag

     angmom_f(1) = angmom_before(1) - angmom_frag(1)
     angmom_f(2) = angmom_before(2) - angmom_frag(2)
     angmom_f(3) = angmom_before(3) - angmom_frag(3)

     mx_frag = 0.0_DP
     my_frag = 0.0_DP
     mz_frag = 0.0_DP
     mvx_frag = 0.0_DP
     mvy_frag = 0.0_DP
     mvz_frag = 0.0_DP
     angmom_frag(1) = 0.0_DP
     angmom_frag(2) = 0.0_DP
     angmom_frag(3) = 0.0_DP
     DO i=1, frags_added
          p_frag(:,i) = p_frag(:,i) + x_f(:)
          vel_frag(:,i) = vel_frag(:,i) + v_f(:)
          veclinmom_after(:) = (m_frag(i) * vel_frag(:,i)) + veclinmom_after(:)
          mx_frag = (p_frag(1,i) * m_frag(i)) + mx_frag
          my_frag = (p_frag(2,i) * m_frag(i)) + my_frag
          mz_frag = (p_frag(3,i) * m_frag(i)) + mz_frag

          mvx_frag = (vel_frag(1,i) * m_frag(i)) + mvx_frag
          mvy_frag = (vel_frag(2,i) * m_frag(i)) + mvy_frag
          mvz_frag = (vel_frag(3,i) * m_frag(i)) + mvz_frag
          call util_crossproduct(p_frag(:,i), vel_frag(:,i), angmom_fragi(:)) 
          angmom_frag(:) = angmom_frag(:) + angmom_fragi(:)
     END DO 

     angmom_com_frag(1) = angmom_before(1) - angmom_frag(1)
     angmom_com_frag(2) = angmom_before(2) - angmom_frag(2)
     angmom_com_frag(3) = angmom_before(3) - angmom_frag(3)

     x_com_frag = mx_frag / SUM(m_frag(:))
     y_com_frag = my_frag / SUM(m_frag(:))
     z_com_frag = mz_frag / SUM(m_frag(:))

     vx_com_frag = mvx_frag / SUM(m_frag(:))
     vy_com_frag = mvy_frag / SUM(m_frag(:))
     vz_com_frag = mvz_frag / SUM(m_frag(:))
     angmom_after = NORM2(angmom_com_frag)
     DL = (angmom_after - NORM2(angmom_before))/ NORM2(angmom_before)
     WRITE(*,*) "util_mom DL/L = ", DL 
     !WRITE(*,*) "util_mom l(1) :", l(1)
     !WRITE(*,*) "util_mom p(1) :", p(1)
     WRITE(*,*) "util_mom p_frag_check :", p_frag_check
     WRITE(*,*) "util_mom v_frag_check :", v_frag_check

     !WRITE(*,*) "util_mom linmom_after: ", NORM2(veclinmom_after)
     WRITE(*,*) "util_mom linmom_diff: ", (NORM2(veclinmom_after) - linmom_before) / linmom_before
     WRITE(*,*) "util_mom x position com diff", (x_com - x_com_frag)
     WRITE(*,*) "util_mom y position com diff", (y_com - y_com_frag)
     WRITE(*,*) "util_mom z position com diff", (z_com - z_com_frag)
     WRITE(*,*) "util_mom x velocity com diff", (vx_com - vx_com_frag)
     WRITE(*,*) "util_mom y velocity com diff", (vy_com - vy_com_frag)
     WRITE(*,*) "util_mom z velocity com diff", (vz_com - vz_com_frag)

   RETURN

END SUBROUTINE util_mom


!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
