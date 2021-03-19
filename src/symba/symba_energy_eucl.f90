subroutine symba_energy_eucl(npl, swiftest_plA, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, te, Ltot, msys)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector and kinetic, potential and total system energy
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
   use swiftest
   use module_interfaces, except_this_one => symba_energy_eucl
   implicit none

! arguments
   integer(I4B), intent(in)                 :: npl
   real(DP), intent(in)                     :: j2rp2, j4rp4
   integer(I4B), dimension(:,:), intent(in) :: k_plpl
   integer(I8B), intent(in)                 :: num_plpl_comparisons
   real(QP), intent(out)                    :: ke, pe, te, msys
   real(QP), dimension(:), intent(out)      :: Ltot
   type(swiftest_pl), intent(inout)         :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   real(QP)                  :: mass, v2, Mcb, rad, Ip, rot2
   real(QP), dimension(NDIM) :: h, x, v, xbcb, rot
   real(DP), dimension(npl)  :: irh
   real(DP)                  :: mtmp, rmag, oblpot
   integer(I8B)              :: k

! executable code

   call coord_h2b(npl, swiftest_plA, mtmp)
   ! Re-do this at quad precision
   msys = swiftest_plA%dMcb + sum(swiftest_plA%mass(2:npl)) + swiftest_plA%Mcb_initial
   ke = 0.0_QP

   x(:) = swiftest_plA%xb(:, 1)
   v(:) = swiftest_plA%vb(:, 1)
   rot(:) = swiftest_plA%rot(:, 1)
   Ip = swiftest_plA%Ip(3, 1)
   mass = swiftest_plA%mass(1)
   rad = swiftest_plA%radius(1)
   v2 = dot_product(v(:), v(:))
   rot2 = dot_product(rot(:), rot(:))
   h(1) = x(2) * v(3) - x(3) * v(2)
   h(2) = x(3) * v(1) - x(1) * v(3)
   h(3) = x(1) * v(2) - x(2) * v(1)
   Ltot(:) = swiftest_plA%Lcb_initial + swiftest_plA%dLcb + mass * h(:)
   ke = 0.5_QP * (swiftest_plA%dMcb + swiftest_plA%Mcb_initial) * (v2 + Ip * rad**2 * rot2)

   !$omp simd private(x,v,v2,mass,h,rot,Ip,rad,rot2) reduction(+:ke,Ltot)
   do i = 2, npl 
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

      Ltot(:) = Ltot(:) + mass * (h(:) + Ip * rot(:) * rad**2)

      ! Kinetic energy from orbit and spin
      ke = ke + 0.5_QP * mass * (v2 + Ip * rad**2 * rot2)
   end do

   ! Do the central body potential energy component first
   xbcb(:) = swiftest_plA%xb(:, 1) 
   Mcb = swiftest_plA%mass(1) 
   pe = 0.0_QP
   !$omp simd private(rmag) reduction(-:pe)
   do i = 2, npl
      rmag = norm2(swiftest_plA%xb(:, i) - xbcb(:))
      pe = pe - Mcb * swiftest_plA%mass(i) / rmag 
   end do

   ! Do the potential energy between pairs of massive bodies
   !$omp parallel do default(private) schedule(static) &
   !$omp shared(num_plpl_comparisons, k_plpl, swiftest_plA) &
   !$omp reduction(-:pe)
   do k = 1, num_plpl_comparisons
       i = k_plpl(1, k)
       j = k_plpl(2, k)
       rmag = norm2(swiftest_plA%xb(:, j) - swiftest_plA%xb(:, i)) 
       if (rmag > tiny(rmag)) pe = pe - swiftest_plA%mass(i) * swiftest_plA%mass(j) / rmag 
   end do
   !$omp end parallel do

   ! Potential energy from the oblateness term
   if (j2rp2 /= 0.0_DP) then
      !$omp simd private(rmag)
      do i = 2, npl
         rmag = norm2(swiftest_plA%xh(:,i))
         irh(i) = 1.0_DP / rmag
      end do
      call obl_pot(npl, swiftest_plA, j2rp2, j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
      pe = pe + oblpot
   end if

   te = ke + pe

   return

end subroutine symba_energy_eucl