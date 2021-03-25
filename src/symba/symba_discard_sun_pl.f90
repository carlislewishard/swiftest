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
   logical             :: ldiscard_this, ldiscard_cb, lescape
   real(DP), dimension(NDIM) :: Lorb0, Lspin0, Lorb1, Lspin1, h, x, v
   real(DP)             :: v2, rot2, Ip, ke0, pe0, rote0, ke1, pe1, rote1

! executable code
   rmin2 = rmin*rmin
   rmax2 = rmax*rmax
   rmaxu2 = rmaxu*rmaxu
   ldiscard_cb = .false.
   lescape = .false.

   !******************* TESTING *************************
   ke0 = 0.0_DP
   pe0 = 0.0_DP
   rote0 = 0.0_DP
   !do i = 1, npl
   !   x(:) = swiftest_plA%xb(:, i)
   !   v(:) = swiftest_plA%vb(:, i)
   !   rot(:) = swiftest_plA%rot(:, i)
   !   Ip = swiftest_plA%Ip(3, i)
   !   mass = swiftest_plA%mass(i)
   !   rad = swiftest_plA%radius(i)
   !   v2 = dot_product(v(:), v(:))
   !   rot2 = dot_product(rot(:), rot(:))
   !   ke = ke + 0.5_DP * mass * v2 
   !   rote = rote + 0.5_DP * mass * Ip * rad**2 * rot2
   !end do
   !!****************************************************

   do i = 2, npl
      if (swiftest_plA%status(i) == ACTIVE) then
         ldiscard_this = .false.
         rh2 = dot_product(swiftest_plA%xh(:,i), swiftest_plA%xh(:,i))
         if ((rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
            ldiscard = .true.
            ldiscard_this = .true.
            lescape = .true.
            swiftest_plA%status(i) = DISCARDED_RMAX
            write(*, *) "Particle ",  swiftest_plA%name(i), " too far from the central body at t = ", t
         else if ((rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
            ldiscard = .true.
            ldiscard_this = .true.
            lescape = .false.
            swiftest_plA%status(i) = DISCARDED_RMIN
            write(*, *) "Particle ", swiftest_plA%name(i), " too close to the central body at t = ", t
         else if (rmaxu >= 0.0_DP) then
            rb2 = dot_product(swiftest_plA%xb(:,i), swiftest_plA%xb(:,i))
            vb2 = dot_product(swiftest_plA%vb(:,i), swiftest_plA%vb(:,i))
            energy = 0.5_DP*vb2 - msys/sqrt(rb2)
            if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
               ldiscard = .true.
               ldiscard_this = .true.
               lescape = .true.
               swiftest_plA%status(i) = DISCARDED_RMAXU
               write(*, *) "Particle ", swiftest_plA%name(i), " is unbound and too far from barycenter at t = ", t
            end if
         end if
         if (ldiscard_this) then
            ldiscard_cb = .true.
            call symba_discard_conserve_mtm(swiftest_plA, i, lescape)
         end if
      end if
   end do

   if (ldiscard_cb) then
      ! Because the central body has changed position, we need to adjust the heliocentric position and velocities of everything
      call coord_b2h(npl, swiftest_plA)
      if (ntp > 0) call coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)
   end if


   !******************* TESTING *************************
   !if (ldiscard_cb) then
   !   Lorb1 = 0.0_DP
   !   Lspin1 = 0.0_DP
   !   do i = 1, npl
   !      if (swiftest_plA%status(i) == DISCARDED_RMIN) cycle
   !      x(:) = swiftest_plA%xb(:, i)
   !      v(:) = swiftest_plA%vb(:, i)
   !      rot(:) = swiftest_plA%rot(:, i)
   !      Ip = swiftest_plA%Ip(3, i)
   !      mass = swiftest_plA%mass(i)
   !      rad = swiftest_plA%radius(i)
   !      v2 = dot_product(v(:), v(:))
   !      rot2 = dot_product(rot(:), rot(:))
   !      h(1) = x(2) * v(3) - x(3) * v(2)
   !      h(2) = x(3) * v(1) - x(1) * v(3)
   !      h(3) = x(1) * v(2) - x(2) * v(1)
   !      Lorb1(:) = Lorb1(:) + mass * h(:)
   !      Lspin1(:) = Lspin1(:) + mass * rad**2 * Ip * rot(:)
   !   end do
   !   write(*,*) 'Before loss: '
   !   write(*,*) 'Lorb   : ', Lorb0(:)
   !   write(*,*) '|Lorb| : ', norm2(Lorb0(:))
   !   write(*,*) 'Lspin  : ', Lspin0(:)
   !   write(*,*) '|Lspin|: ', norm2(Lspin0(:))
   !   write(*,*) 'Ltot   : ', Lspin0(:) + Lorb0(:)
   !   write(*,*) '|Ltot| : ', norm2(Lspin0(:) + Lorb0(:))
!
!      write(*,*) 'After loss: '
!      write(*,*) 'Lorb   : ', Lorb1(:)
!      write(*,*) '|Lorb| : ', norm2(Lorb1(:))
!      write(*,*) 'Lspin  : ', Lspin1(:)
!      write(*,*) '|Lspin|: ', norm2(Lspin1(:))
!      write(*,*) 'Ltot   : ', Lspin1(:) + Lorb1(:)
!      write(*,*) '|Ltot| : ', norm2(Lspin1(:) + Lorb1(:))
!
!      write(*,*) 'Difference (after - before): '
!      write(*,*) 'Lorb1 - Lorb0      : ', Lorb1(:) - Lorb0(:)
!      write(*,*) '|Lorb1| - |Lorb0|  : ', norm2(Lorb1(:)) - norm2(Lorb0(:))
!      write(*,*) 'Lspin1 - Lspin0    : ', Lspin1(:) - Lspin0(:)
!      write(*,*) '|Lspin1| - |Lspin0|: ', norm2(Lspin1(:)) - norm2(Lspin0(:))
!      write(*,*) 'Ltot1 - Ltot0      : ', Lspin1(:) + Lorb1(:) - Lspin0(:) - Lorb0(:)
!      write(*,*) '|Ltot1| - |Ltot0|  : ', norm2(Lspin1(:) + Lorb1(:)) - norm2(Lspin0(:) + Lorb0(:))
!
!      write(*,*) 'Ratio (before / after): '
!      write(*,*) 'Lorb0 / Lorb1      : ', Lorb0(:) / Lorb1(:)
!      write(*,*) '|Lorb0| / |Lorb1|  : ', norm2(Lorb0(:)) / norm2(Lorb1(:))
!      write(*,*) 'Lspin0 / Lspin1    : ', Lspin0(:) / Lspin1(:)
!      write(*,*) '|Lspin0| / |Lspin1|: ', norm2(Lspin0(:)) / norm2(Lspin1(:))
!      write(*,*) 'Ltot0 / Ltot1      : ', (Lspin0(:) + Lorb0(:)) / (Lspin1(:) + Lorb1(:))
!      write(*,*) '|Ltot0| / |Ltot1|  : ', norm2(Lspin0(:) + Lorb0(:)) / norm2(Lspin1(:) + Lorb1(:))
!      read(*,*)
!   end if
   !****************************************************


   return

end subroutine symba_discard_sun_pl