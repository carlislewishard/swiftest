subroutine symba_energy(npl, swiftest_plA, j2rp2, j4rp4, ke, pe, te, htot)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector and kinetic, potential and total system energy
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
   real(DP), intent(out)            :: ke, pe, te
   real(DP), dimension(:), intent(out) :: htot
   type(swiftest_pl), intent(inout)     :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   real(DP)                  :: mass, msys, rmag, v2, oblpot
   real(DP), dimension(NDIM) :: h, dx, x, v
   real(DP), dimension(npl)  :: irh


! executable code

   call coord_h2b(npl, swiftest_plA, msys)
   htot = 0.0_DP
   ke = 0.0_DP


   !$omp simd private(x,v,v2,mass,h) reduction(+:ke,htot)
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

   pe = 0.0_DP
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
         irh(i) = 1.0_DP / rmag
      end do
      call obl_pot(npl, swiftest_plA, j2rp2, j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
      pe = pe + oblpot
   end if

   te = ke + pe

   return

end subroutine symba_energy