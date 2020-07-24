!**********************************************************************************************************************************
!
!  Unit Name   : drift_one
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Perform Danby drift for one body, redoing drift with smaller substeps if original accuracy is insufficient
!
!  Input
!    Arguments : mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                x     : position of body to drift
!                v     : velocity of body to drift
!                dt    : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : x     : position of body to drift
!                v     : velocity of body to drift
!                iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_one(mu, x, v, dt, iflag)
!
!  Notes       : Adapted from Hal Levison and Martin Duncan's Swift routine drift_one.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_one(mu, x, v, dt, iflag, n)

! Modules
     USE swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => drift_one
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), dimension(:), INTENT(OUT) :: iflag
     REAL(DP), INTENT(IN)                     :: mu, dt
     REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: x, v
     integer(I4B), intent(in)                :: n

! Internals
     INTEGER(I4B) :: i, k
     real(DP), dimension(:), allocatable :: muvec, dttmp
     real(DP), dimension(:), allocatable :: px, py, pz, vx, vy, vz
     logical(DP), dimension(:), allocatable :: badrun

! Executable code
     allocate(muvec(n))
     allocate(dttmp(n))
     allocate(px(n))
     allocate(py(n))
     allocate(pz(n))
     allocate(vx(n))
     allocate(vy(n))
     allocate(vz(n))
     allocate(badrun(n))
     muvec(:) = mu
     dttmp(:) = dt
     px(:) = x(1, :)
     py(:) = x(2, :)
     pz(:) = x(3, :)
     vx(:) = v(1, :)
     vy(:) = v(2, :)
     vz(:) = v(3, :)
     do i = 1, n
         CALL drift_dan(muvec(i), px(i), py(i), pz(i), vx(i), vy(i), vz(i), dttmp(i), iflag(i))
     end do
     badrun(:) = (iflag(:) /= 0)
     if (any(badrun(:))) then
         muvec(:) = pack(muvec(:), badrun(:))
         px(:) = pack(px(:), badrun(:))
         py(:) = pack(py(:), badrun(:))
         pz(:) = pack(pz(:), badrun(:))
         vx(:) = pack(vx(:), badrun(:))
         vy(:) = pack(vy(:), badrun(:))
         vz(:) = pack(vz(:), badrun(:))
         dttmp(:) = 0.1_DP * pack(dttmp(:), badrun(:))
         iflag(:) = pack(iflag(:), badrun(:))
         do k = 1, 10
            do i = 1, count(badrun(:))
               call drift_dan(muvec(i), px(i), py(i), pz(i), vx(i), vy(i), vz(i), dttmp(i), iflag(i))
            end do
         end do
         if (all(iflag(:) == 0)) then
            x(1, :) = unpack(px(:), badrun(:), x(1, :)) 
            x(2, :) = unpack(py(:), badrun(:), x(2, :)) 
            x(3, :) = unpack(pz(:), badrun(:), x(3, :)) 
            v(1, :) = unpack(vx(:), badrun(:), v(1, :)) 
            v(2, :) = unpack(vy(:), badrun(:), v(2, :)) 
            v(3, :) = unpack(vz(:), badrun(:), v(3, :)) 
            return
         end if
     end if
      x(1, :) = px(:)
      x(2, :) = py(:)
      x(3, :) = pz(:)
      v(1, :) = vx(:)
      v(2, :) = vy(:)
      v(3, :) = vz(:)

     RETURN

END SUBROUTINE drift_one
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
