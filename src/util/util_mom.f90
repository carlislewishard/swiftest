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
SUBROUTINE util_mom(m1, xh1, vb1, m2, xh2, vb2, frags_added, nstart, m_frag, r_circle, theta, p_frag, vel_frag)

!x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag)

! Modules
   USE swiftest_globals
   USE swiftest_data_structures
   USE module_interfaces, EXCEPT_THIS_ONE => util_mom
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                               :: frags_added, nstart
   REAL(DP), INTENT(IN)                                   :: m1, m2, r_circle, theta
   REAL(DP), DIMENSION(:), INTENT(IN)                     :: xh1, vb1, xh2, vb2
   REAL(DP), DIMENSION(:), INTENT(IN)                     :: m_frag
   REAL(DP), DIMENSION(:, :), INTENT(OUT)                 :: p_frag, vel_frag

! Internals

   INTEGER(I4B)                                           :: i 
   REAL(DP)                                               :: v_col, A, v2el, v2esc
   REAL(DP)                                               :: linmom_before,  linmom_after, DL
   REAL(DP), DIMENSION(NDIM)                              :: veclinmom_after, xhvb1, xhvb2, vecangmom_after, xvrel
   REAL(DP)                                               :: p_frag_check, v_frag_check, B, m_frag_tot
   REAL(DP), DIMENSION(NDIM)                              :: xr, l, kk, p, v_f, x_f, angmom_frag, angmom_fragi, angmom_com_frag
   REAL(DP), DIMENSION(NDIM)                              :: angmom_f, angmom_before, angmom_after
   integer(I4B), save                                     :: thetashift = 0
   integer(I4B), parameter                                :: SHIFTMAX = 9
   real(DP)                                               :: phase_ang
   real(DP), dimension(NDIM)                              :: p_com, v_com, mp_frag, mv_frag, p_com_frag, v_com_frag

! Executable code

   ! Shifts the starting circle of fragments around so that multiple fragments generated in from a single body in a single time step 
   ! don't pile up on top of each other
   phase_ang = theta * thetashift / SHIFTMAX
   thetashift = thetashift + 1
   if (thetashift >= shiftmax) thetashift = 0


   linmom_before = NORM2(m1*vb1(:) + m2*vb2(:))
   call util_crossproduct(xh1,vb1,xhvb1)
   call util_crossproduct(xh2,vb2,xhvb2)
   angmom_before = (m1*xhvb1+m2*xhvb2)


        ! Find COM
   p_com(:) = ((xh1(:) * m1) + (xh2(:) * m2)) / (m1 + m2)
   v_com(:) = ((vb1(:) * m1) + (vb2(:) * m2)) / (m1 + m2)

   ! Find Collision velocity
     v_col = NORM2(vb2(:) - vb1(:))
     xr(:) = xh2(:) - xh1(:)
     l(:) = (vb2(:) - vb1(:)) 
     call util_crossproduct(l,xr,xvrel)
     kk(:) = xvrel !angmom_before
     !call util_crossproduct(l(:), xr(:), kk(:))
     call util_crossproduct(kk(:), l(:), p(:))
     kk(:) = kk(:) / NORM2(kk(:))
     l(:) = l(:) / NORM2(l(:))
     p(:) = p(:) / NORM2(p(:)) 

     veclinmom_after(:) = 0.0_DP
     vecangmom_after(:) = 0.0_DP

     mp_frag = 0.0_DP
     mv_frag = 0.0_DP
     p_frag_check = 0.0_DP
     v_frag_check = 0.0_DP
     angmom_frag(:) = 0.0_DP
     v2esc = 2.0_DP * GC * (m1+m2) / (NORM2(xr))
     v2el = v2esc - 2.0_DP*(m1+m2)*GC*(1.0_DP/(NORM2(xr)) - 1.0_DP/r_circle)
     A = - (SQRT(v2el))

     !WRITE(*,*) "UTIL_MOM A", A

     B = r_circle
     DO i=1, frags_added
          p_frag(:,i) = (- B  * cos(phase_ang + theta * i))*l(:) + (- B  * sin(phase_ang + theta * i))*p(:) + p_com(:)
          p_frag_check = - (B * cos(phase_ang + theta * i)) + p_frag_check
          vel_frag(:,i) = (((A * cos(phase_ang + theta * i))*l(:)) + ((A * sin(phase_ang + theta * i))*p(:)))  + v_com(:)

          v_frag_check = (A * cos(phase_ang + theta * i)) + v_frag_check

          mp_frag = (p_frag(:,i) * m_frag(i)) + mp_frag(:)

          mv_frag = (vel_frag(:,i) * m_frag(i)) + mv_frag(:)

          call util_crossproduct(p_frag(:,i), vel_frag(:,i), angmom_fragi(:)) 
          angmom_frag(:) = angmom_frag(:) + angmom_fragi(:)
     END DO

     m_frag_tot = SUM(m_frag(:))
     p_com_frag(:) = mp_frag(:) / m_frag_tot
     v_com_frag(:) = mv_frag(:) / m_frag_tot
     x_f(:) =  p_com(:) - p_com_frag(:)
     v_f(:) = v_com(:) - v_com_frag(:)
     angmom_f(:) = angmom_before(:) - angmom_frag(:)

     mp_frag(:) = 0.0_DP
     mv_frag(:) = 0.0_DP
     angmom_frag(:) = 0.0_DP

     DO i=1, frags_added
          p_frag(:,i) = p_frag(:,i) + x_f(:)
          vel_frag(:,i) = vel_frag(:,i) + v_f(:)
          veclinmom_after(:) = (m_frag(i) * vel_frag(:,i)) + veclinmom_after(:)
          mp_frag(:) = (p_frag(:,i) * m_frag(i)) + mp_frag(:)
          mv_frag(:) = (vel_frag(:,i) * m_frag(i)) + mv_frag(:)
          call util_crossproduct(p_frag(:,i)*m_frag(i), vel_frag(:,i), angmom_fragi(:)) 
          angmom_frag(:) = angmom_frag(:) + angmom_fragi(:)
     END DO 

     p_com_frag(:) = mp_frag(:) / SUM(m_frag(:))
     v_com_frag(:) = mv_frag(:) / SUM(m_frag(:))
     angmom_after = angmom_frag(:)
     DL = NORM2(angmom_after(:) - angmom_before(:)) / NORM2(angmom_before(:))

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
