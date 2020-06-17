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
SUBROUTINE util_mom(m1, x1, v1, m2, x2, v2, frags_added, nstart, m_frag, r_circle, theta, p_frag, vel_frag)

!x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag)

! Modules
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_interfaces, EXCEPT_THIS_ONE => util_mom
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: frags_added, nstart
   REAL(DP), INTENT(IN)                                   :: m1, m2, r_circle, theta
   REAL(DP), DIMENSION(NDIM), INTENT(IN)                  :: x1, v1, x2, v2
   REAL(DP), DIMENSION(frags_added), INTENT(IN)           :: m_frag
   REAL(DP), DIMENSION(NDIM, frags_added), INTENT(INOUT)  :: p_frag, vel_frag

! Internals

   INTEGER(I4B)                                           :: i 
   REAL(DP)                                               :: x_com, y_com, z_com, vx_com, vy_com, vz_com, v_col, A
   REAL(DP)                                               :: linmom_before
   REAL(DP), DIMENSION(NDIM)                              :: linmom_after
   REAL(DP)                                               :: mx_frag, my_frag, mz_frag, mvx_frag, mvy_frag, mvz_frag
   REAL(DP)                                               :: x_com_frag, y_com_frag, z_com_frag, vx_com_frag, vy_com_frag
   REAL(DP)                                               :: vz_com_frag, p_frag_check, v_frag_check
   REAL(DP), DIMENSION(NDIM)                              :: xr, l, kk, p

! Executable code

     linmom_before = NORM2(m1*v1(:) + m2*v2(:))

     WRITE(*,*) "util_mom linmom_before: ", linmom_before

        ! Find COM
   x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! Find Collision velocity
      v_col = NORM2(v2(:) - v1(:))

     xr(:) = x2(:) - x1(:)
     l(:) = (v2(:) - v1(:)) / NORM2(v2(:)-v1(:))
     call util_crossproduct(xr(:), l(:), p(:))
     p(:) = p(:) / NORM2(p(:))
     call util_crossproduct(l(:), p(:), kk(:))
     kk(:) = kk(:) / NORM2(kk(:)) 

     linmom_after = 0.0_DP

     WRITE(*,*) "util_mom check 1"

     mx_frag = 0.0_DP
     my_frag = 0.0_DP
     mz_frag = 0.0_DP
     mvx_frag = 0.0_DP
     mvy_frag = 0.0_DP
     mvz_frag = 0.0_DP
     p_frag_check = 0.0_DP
     v_frag_check = 0.0_DP

     DO i=1, frags_added
         WRITE(*,*) "util_mom check 2"
         WRITE(*,*) "m_frag(i)", m_frag(i)

          A = (v_col * (m1 + m2))

          WRITE(*,*) "util_mom check 3"
          
          p_frag(1,i) = (r_circle / m_frag(i) * cos(theta * i))*l(1) + (r_circle / m_frag(i) * sin(theta * i))*p(1) + x_com
          p_frag(2,i) = (r_circle / m_frag(i) * cos(theta * i))*l(2) + (r_circle / m_frag(i) * sin(theta * i))*p(2) + y_com
          p_frag(3,i) = (r_circle / m_frag(i) * cos(theta * i))*l(3) + (r_circle / m_frag(i) * sin(theta * i))*p(3) + z_com

          p_frag_check = (r_circle * cos(theta * i)) + p_frag_check

          WRITE(*,*) "util_mom check 4"

          vel_frag(1,i) = (((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1))) / (m_frag(i)) + vx_com
          vel_frag(2,i) = (((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2))) / (m_frag(i)) + vy_com
          vel_frag(3,i) = (((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3))) / (m_frag(i)) + vz_com

          v_frag_check = (A * cos(theta * i)) + v_frag_check

          linmom_after(:) = (m_frag(i) * vel_frag(:,i)) + linmom_after

          mx_frag = (p_frag(1,i) * m_frag(i)) + mx_frag
          my_frag = (p_frag(2,i) * m_frag(i)) + my_frag
          mz_frag = (p_frag(3,i) * m_frag(i)) + mz_frag

          mvx_frag = (vel_frag(1,i) * m_frag(i)) + mvx_frag
          mvy_frag = (vel_frag(2,i) * m_frag(i)) + mvy_frag
          mvz_frag = (vel_frag(3,i) * m_frag(i)) + mvz_frag

     END DO
 
     WRITE(*,*) "util_mom l(1) :", l(1)
     WRITE(*,*) "util_mom p(1) :", p(1)
     WRITE(*,*) "util_mom p_frag_check :", p_frag_check
     WRITE(*,*) "util_mom v_frag_check :", v_frag_check

     WRITE(*,*) "util_mom check 5"

     x_com_frag = mx_frag / SUM(m_frag(:))
     y_com_frag = my_frag / SUM(m_frag(:))
     z_com_frag = mz_frag / SUM(m_frag(:))

     vx_com_frag = mvx_frag / SUM(m_frag(:))
     vy_com_frag = mvy_frag / SUM(m_frag(:))
     vz_com_frag = mvz_frag / SUM(m_frag(:))

     WRITE(*,*) "util_mom linmom_after: ", NORM2(linmom_after)
     WRITE(*,*) "util_mom linmom_diff: ", (NORM2(linmom_after) - linmom_before) / linmom_before
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
