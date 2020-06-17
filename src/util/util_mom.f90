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
SUBROUTINE util_mom(m1, x1, v1, m2, x2, v2, frags_added, nstart, m_frag, r_circle, theta, x_frag, v_frag)

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
   REAL(DP), DIMENSION(NDIM, frags_added), INTENT(INOUT)  :: x_frag, v_frag

! Internals

   INTEGER(I4B)                                           :: i 
   REAL(DP)                                               :: x_com, y_com, z_com, vx_com, vy_com, vz_com, v_col, A
   REAL(DP), DIMENSION(NDIM)                              :: xr, l, kk, p

! Executable code

   !TEMPORARY
   
   !interface 
   !   function cross_product_mom(ar1,ar2) result(ans)
   !      use swiftest
   !      implicit none
   !      real(DP),dimension(3),intent(in) :: ar1,ar2
   !      real(DP),dimension(3)             :: ans
   !   end function cross_product_mom
   !end interface

   ! Find COM
   !x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   !y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   !z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   !vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   !vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   !vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! Find Collision velocity
   v_col = NORM2(v2(:) - v1(:))

   DO i=1, frags_added

     xr(:) = x2(:) - x1(:)
     l(:) = (v2(:) - v1(:)) / NORM2(v2(:)-v1(:))
     p(:) = cross_product_mom(xr(:) / NORM2(xr(:)), l(:))
     kk(:) = cross_product_mom(l(:),p(:))

     DO i=1, frags_added

          A = v_col * (m1 + m2) * (1.0_DP / mergeadd_list%mass(nstart + i))
          
          x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
          y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
          z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

          vx_frag = ((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1)) + vx_com
          vy_frag = ((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2)) + vy_com
          vz_frag = ((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3)) + vz_com
     END DO

   RETURN

END SUBROUTINE util_mom


!function cross_product_mom(ar1,ar2) result(ans)
!   use swiftest
!   implicit none
!   
!   real(DP),dimension(3),intent(in) :: ar1,ar2
!   real(DP),dimension(3)             :: ans

 !  ans(1) = ar1(2) * ar2(3) - ar1(3) * ar2(2)
 !  ans(2) = ar1(3) * ar2(1) - ar1(1) * ar2(3)
 !  ans(3) = ar1(1) * ar2(2) - ar1(2) * ar2(1)

  ! return 
!end function cross_product_mom

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
