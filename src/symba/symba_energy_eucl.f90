subroutine symba_energy_eucl(npl, swiftest_plA, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, te, htot)
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
   integer(I4B), intent(in)         :: npl
   real(DP), intent(in)             :: j2rp2, j4rp4
   integer(I4B), dimension(:,:), intent(in) :: k_plpl
   integer(I8B), intent(in)         :: num_plpl_comparisons
   real(DP), intent(out)            :: ke, pe, te
   real(DP), dimension(:), intent(out) :: htot
   type(swiftest_pl), intent(inout)     :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   real(DP)                  :: mass, msys, rmag, v2, oblpot
   real(DP), dimension(NDIM) :: h, dx, x, v
   real(DP), dimension(npl)  :: irh
   integer(I8B)              :: k

! executable code

   call coord_h2b(npl, swiftest_plA, msys)
   htot = 0.0_DP
   ke = 0.0_DP

   ! Kinetic energy of all bodies, inclusing the central body
   !$omp simd
   do i = 1, npl 
      x(:) = swiftest_plA%xb(:, i)
      v(:) = swiftest_plA%vb(:, i)
      mass = swiftest_plA%mass(i)
      h(1) = mass * (x(2) * v(3) - x(3) * v(2))
      h(2) = mass * (x(3) * v(1) - x(1) * v(3))
      h(3) = mass * (x(1) * v(2) - x(2) * v(1))
      htot(:) = htot(:) + h(:)
      v2 = dot_product(v(:), v(:))
      ke = ke + 0.5_DP * mass * v2
   end do

   ! Do the central body potential energy component first
   pe = 0.0_DP
   !$omp simd
   do i = 2, npl
      rmag = norm2(swiftest_plA%xb(:, i) - swiftest_plA%xb(:, 1))
      pe = pe - swiftest_plA%mass(i) * swiftest_plA%mass(1) / rmag 
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
      !$omp simd
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