!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step the fluid ring in time
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                t              : time
!                npl            : number of planets
!                nplmax         : maximum allowed number of planets
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                eoffset        : energy offset (net energy lost to ring)
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_step(lfirst, t, npl, nplmax, ntp, ntpmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt)
!
!  Notes       : Adapted from Andy Hesselbrock's RING-MOONS Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_step(lfirst, t, rmin,npl, nplmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt,ring)

! Modules
     USE module_parameters
     USE module_symba
     USE module_ringmoons
     USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
     INTEGER(I4B), INTENT(IN)                         :: npl, nplmax
     REAL(DP), INTENT(IN)                             :: t, rmin,j2rp2, j4rp4, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     TYPE(symba_pl), POINTER                          :: symba_pl1P
     TYPE(ringmoons_ring),INTENT(INOUT) :: ring
! Internals
      integer(I4B) :: i

! Executable code
      call ringmoons_pde_solver(symba_pl1P%helio%swifter%mass,rmin,dt,ring)
      
      open(unit=22,file='test.out',status='replace')
      do i = 1,ring%N
         write(22,*) ring%r(i) * DU2CM,ring%nu(i) * DU2CM**2/TU2S,ring%sigma(i)*MU2GM/DU2CM**2
      end do
      close(22)
      read(*,*)

     RETURN

END SUBROUTINE ringmoons_step
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
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
