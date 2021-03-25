!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_sun_pl
!  Unit Type   : subroutine
!  Project   : Swiftest
!  Package   : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their positions relative to the Sun
!
!  Input
!    Arguments : t        : time
!          npl      : number of planets
!          msys       : total system mass
!          swifter_pl1P : pointer to head of Swifter planet structure linked-list
!          rmin       : minimum allowed heliocentric radius
!          rmax       : maximum allowed heliocentric radius
!          rmaxu      : maximum allowed heliocentric radius for unbound planets
!          ldiscard    : logical flag indicating whether any planets are discarded
!    Terminal  : none
!    File    : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!          ldiscard    : logical flag indicating whether any planets are discarded
!    Terminal  : status message
!    File    : none
!
!  Invocation  : CALL symba_discard_sun_pl(t, npl, msys, swifter_pl1P, rmin, rmax, rmaxu, ldiscard)
!
!  Notes     : Adapted from Hal Levison's Swift routine discard_massive5.f
!
!**********************************************************************************************************************************
subroutine symba_discard_sun_pl(t, npl, ntp, msys, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, ldiscard)
   !! author: David A. Minton
   !!
   !! Check to see if planets should be discarded based on their positions relative to the Sun.
   !! Updates the mass and angular momentum of the central body
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f

! modules
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => symba_discard_sun_pl
   implicit none

! arguments
   integer(I4B), intent(in)       :: npl, ntp
   real(DP), intent(in)         :: t, msys, rmin, rmax, rmaxu
   type(swiftest_pl), intent(inout) :: swiftest_plA
   type(swiftest_tp), intent(inout) :: swiftest_tpA
   logical(LGT), intent(inout)    :: ldiscard

! internals
   integer(I4B)          :: i
   real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
   real(DP)            :: mass, rad, Ipz, Ipcbz, Mcb, radcb
   logical             :: ldiscard_this, ldiscard_cb
   real(DP), dimension(NDIM) :: Lpl, Lcb, rot, rotcb, xb, vb, xbcb, vbcb, xcom, vcom

! executable code
   rmin2 = rmin*rmin
   rmax2 = rmax*rmax
   rmaxu2 = rmaxu*rmaxu
   ldiscard_cb = .false.
   do i = 2, npl
      if (swiftest_plA%status(i) == ACTIVE) then
         ldiscard_this = .false.
         rh2 = dot_product(swiftest_plA%xh(:,i), swiftest_plA%xh(:,i))
         if ((rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
            ldiscard = .true.
            ldiscard_this = .true.
            swiftest_plA%status(i) = DISCARDED_RMAX
            swiftest_plA%Mescape = swiftest_plA%Mescape + mass
            write(*, *) "Particle ",  swiftest_plA%name(i), " too far from the central body at t = ", t
         else if ((rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
            ldiscard = .true.
            ldiscard_this = .true.
            swiftest_plA%status(i) = DISCARDED_RMIN
            write(*, *) "Particle ", swiftest_plA%name(i), " too close to the central body at t = ", t
         else if (rmaxu >= 0.0_DP) then
            rb2 = dot_product(swiftest_plA%xb(:,i), swiftest_plA%xb(:,i))
            vb2 = dot_product(swiftest_plA%vb(:,i), swiftest_plA%vb(:,i))
            energy = 0.5_DP*vb2 - msys/sqrt(rb2)
            if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
               ldiscard = .true.
               ldiscard_this = .true.
               swiftest_plA%status(i) = DISCARDED_RMAXU
               swiftest_plA%Mescape = swiftest_plA%Mescape + mass
               write(*, *) "Particle ", swiftest_plA%name(i), " is unbound and too far from barycenter at t = ", t
            end if
         end if
         if (ldiscard_this) then
            ldiscard_cb = .true.
            xb(:) = swiftest_plA%xb(:,i)
            vb(:) = swiftest_plA%vb(:,i)
            rot(:) = swiftest_plA%rot(:,i)
            Ipz = swiftest_plA%Ip(3,i)
            rad = swiftest_plA%radius(i)
            mass = swiftest_plA%mass(i) 

            xbcb(:) = swiftest_plA%xb(:,1)
            vbcb(:) = swiftest_plA%vb(:,1)
            rotcb(:) = swiftest_plA%rot(:,1)
            Mcb = swiftest_plA%mass(1)
            Ipcbz = swiftest_plA%Ip(3,1)
            radcb = swiftest_plA%radius(1)
            xcom(:) = (mass * xb(:) + Mcb * xbcb(:)) / (mass + Mcb)
            vcom(:) = (mass * vb(:) + Mcb * vbcb(:)) / (mass + Mcb)

            call util_crossproduct(xb-xcom,vb-vcom,Lpl)
            Lpl(:) = mass * (Lpl(:) + rad**2 * Ipz * rot(:))

            call util_crossproduct(xbcb-xcom,vbcb-vcom,Lcb)
            Lcb(:) = Mcb * Lcb(:) 
 
            ! Add planet mass to central body accumulator
            swiftest_plA%dMcb = swiftest_plA%dMcb + mass

            ! Update mass of central body to be consistent with its total mass
            Mcb = swiftest_plA%Mcb_initial + swiftest_plA%dMcb
            swiftest_plA%mass(1) = Mcb
            
            ! Add planet angular momentum to central body accumulator
            swiftest_plA%dLcb(:) = swiftest_plA%dLcb(:) + Lpl(:) + Lcb(:) 

            ! Update rotation of central body to by consistent with its angular momentum 
            swiftest_plA%rot(:,1) = (swiftest_plA%Lcb_initial(:) + swiftest_plA%dLcb(:)) / (Ipcbz * Mcb * radcb**2)        
             
            ! Update position and velocity of central body
            swiftest_plA%xb(:,1) = xcom(:)
            swiftest_plA%vb(:,1) = vcom(:)
         end if
      end if
   end do

   if (ldiscard_cb) then
      ! Because the central body has changed position, we need to adjust the heliocentric position and velocities of everything
      call coord_b2h(npl, swiftest_plA)
      if (ntp > 0) call coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)
   end if


   return

end subroutine symba_discard_sun_pl