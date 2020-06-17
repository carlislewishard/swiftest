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
   REAL(DP)                                               :: x_com, y_com, z_com, vx_com, vy_com, vz_com, v_col, A, linmom_before
   REAL(DP), DIMENSION(NDIM)                              :: xr, l, kk, p

! Executable code

     linmom_before = m1*v1(:) + m2*v2(:)

     WRITE(*,*) "util_mom linmom_before: ", linmom_before

   ! Find Collision velocity
   v_col = NORM2(v2(:) - v1(:))

   DO i=1, frags_added

     xr(:) = x2(:) - x1(:)
     l(:) = (v2(:) - v1(:)) / NORM2(v2(:)-v1(:))
     call util_crossproduct(xr(:), l(:), p(:))
     p(:) = p(:) / NORM2(p(:))
     call util_crossproduct(l(:), p(:), kk(:))
     kk(:) = kk(:) / NORM2(kk(:)) 

     linmom_after = 0.0_DP

     DO i=1, frags_added

          A = v_col * (m1 + m2) * (1.0_DP / mergeadd_list%mass(nstart + i))
          
          p_frag(1,i) = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
          p_frag(2,i) = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
          p_frag(3,i) = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

          vel_frag(1,i) = ((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1)) + vx_com
          vel_frag(2,i) = ((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2)) + vy_com
          vel_frag(3,i) = ((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3)) + vz_com

          linmom_after = (m_frag(i) * vel_frag(:,i)) + linmom_after
     END DO



     WRITE(*,*) "util_mom linmom_after: ", linmom_after
     WRITE(*,*) "util_mom linmom_diff: ", (linmom_after - linmom_before) / linmom_before

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
