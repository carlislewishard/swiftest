subroutine symba_energy(npl, swiftest_plA, j2rp2, j4rp4, ke, pe, te, Ltot, msys)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector, kinetic, potential and total 
   !! system energy, and total system mass
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
   use swiftest
   use module_interfaces, except_this_one => symba_energy
   implicit none

! arguments
   integer(I4B), intent(in)         :: npl
   real(DP), intent(in)             :: j2rp2, j4rp4
   real(QP), intent(out)            :: ke, pe, te, msys
   real(QP), dimension(:), intent(out) :: Ltot
   type(swiftest_pl), intent(inout)     :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   real(QP)                  :: mass, v2, Ip, rot2, rad
   real(QP), dimension(NDIM) :: h, dx, x, v, rot
   real(DP), dimension(npl)  :: irh
   real(DP)                  :: mtmp, rmag, oblpot


! executable code

   call coord_h2b(npl, swiftest_plA, mtmp)
   ! Re-do this at quad precision
   msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl)) + swiftest_plA%Mcb_initial
   Ltot = 0.0_QP
   ke = 0.0_QP

   do i = 1, npl
      x(:) = swiftest_plA%xb(:, i)
      v(:) = swiftest_plA%vb(:, i)
      rot(:) = swiftest_plA%rot(:, i)
      Ip = swiftest_plA%Ip(3, i)
      mass = swiftest_plA%mass(i)
      rad = swiftest_plA%radius(i)
      v2 = dot_product(v(:), v(:))
      rot2 = dot_product(rot(:), rot(:))
      h(1) = x(2) * v(3) - x(3) * v(2)
      h(2) = x(3) * v(1) - x(1) * v(3)
      h(3) = x(1) * v(2) - x(2) * v(1)

      ! Angular momentum from orbit and spin
      if (i ==1) then ! Use higher precision angular momentum tracker for the central body
         Ltot(:) = swiftest_plA%Lcb_initial + swiftest_plA%dLcb + mass * h(:)
      else
         Ltot(:) = Ltot(:) + mass * (h(:) + Ip * rot(:) * rad**2)
      end if

      ! Kinetic energy from orbit and spin
      ke = ke + 0.5_QP * mass * (v2 + Ip * rad**2 * rot2)
   end do

   pe = 0.0_QP
   !$omp parallel do default(private) &
   !$omp shared (swiftest_plA, npl) &
   !$omp reduction (-:pe)
   do i = 1, npl - 1
      do j = i + 1, npl
         dx(:) = swiftest_plA%xb(:, j) - swiftest_plA%xb(:, i) 
         rmag = norm2(dx(:)) 
         if (rmag > tiny(rmag)) pe = pe - swiftest_plA%mass(i) * swiftest_plA%mass(j) / rmag 
      end do
   end do
   !$omp end parallel do

   if (j2rp2 /= 0.0_DP) then
      !$omp simd private(rmag)
      do i = 2, npl
         rmag = norm2(swiftest_plA%xh(:,i))
         irh(i) = 1.0_QP / rmag
      end do
      call obl_pot(npl, swiftest_plA, j2rp2, j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
      pe = pe + oblpot
   end if

   te = ke + pe

   return

end subroutine symba_energy