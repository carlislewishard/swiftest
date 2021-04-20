subroutine symba_energy_eucl(npl, swiftest_plA, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, te, Ltot, msys)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector and kinetic, potential and total system energy
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => symba_energy_eucl
   implicit none

! arguments
   integer(I4B), intent(in)                 :: npl
   real(DP), intent(in)                     :: j2rp2, j4rp4
   integer(I4B), dimension(:,:), intent(in) :: k_plpl
   integer(I8B), intent(in)                 :: num_plpl_comparisons
   real(DP), intent(out)                    :: ke, pe, te, msys
   real(DP), dimension(:), intent(out)      :: Ltot
   type(swiftest_pl), intent(inout)         :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   integer(I8B)              :: k
   real(DP)                  :: rmag, v2, rot2, oblpot
   real(DP), dimension(NDIM) :: h
   real(DP), dimension(npl)  :: irh, kepl
   real(DP), dimension(NDIM, npl) :: Lpl 
   real(DP), dimension(num_plpl_comparisons) :: pepl 

! executable code

   call coord_h2b(npl, swiftest_plA, msys)
   Ltot = 0.0_DP
   ke = 0.0_DP
   associate(xb => swiftest_plA%xb, vb => swiftest_plA%vb, mass => swiftest_plA%mass, radius => swiftest_plA%radius, &
      Ip => swiftest_plA%Ip, rot => swiftest_plA%rot, xh => swiftest_plA%xh, status => swiftest_plA%status)
      kepl(:) = 0.0_DP
      Lpl(:,:) = 0.0_DP
      !$omp simd private(v2, rot2, h)
      do i = 1, npl
         if (status(i) /= ACTIVE) cycle
         v2 = dot_product(vb(:,i), vb(:,i))
         rot2 = dot_product(rot(:,i), rot(:,i))
         h(1) = xb(2,i) * vb(3,i) - xb(3,i) * vb(2,i)
         h(2) = xb(3,i) * vb(1,i) - xb(1,i) * vb(3,i)
         h(3) = xb(1,i) * vb(2,i) - xb(2,i) * vb(1,i)

         ! Angular momentum from orbit and spin
         Lpl(1, i) = Ip(3,i) * radius(i)**2 * rot(1,i) + h(1)
         Lpl(2, i) = Ip(3,i) * radius(i)**2 * rot(2,i) + h(2)
         Lpl(3, i) = Ip(3,i) * radius(i)**2 * rot(3,i) + h(3)
         Lpl(:,i) = mass(i) * Lpl(:,i)

         ! Kinetic energy from orbit and spin
         kepl(i) = 0.5_DP * mass(i) * (Ip(3,i) * radius(i)**2 * rot2 + v2)
      end do

      ke = sum(kepl(:))
      Ltot(:) =  sum(Lpl(:,:), dim=2)

      ! Do the central body potential energy component first
      pe = 0.0_DP
      !$omp simd reduction(-:pe)
      do i = 2, npl
         if (status(i) == ACTIVE) pe = pe - mass(1) * mass(i) / norm2(xh(:, i))
      end do

      ! Do the potential energy between pairs of massive bodies
      pepl(:) = 0.0_DP
      !$omp parallel do default(private) schedule(auto) &
      !$omp shared(num_plpl_comparisons, k_plpl, swiftest_plA, pepl)
      do k = 1, num_plpl_comparisons
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         rmag = norm2(xb(:, j) - xb(:, i)) 
         if (rmag > tiny(rmag) .and. (status(i) == ACTIVE) .and. (status(j) == ACTIVE)) then
            pepl(k) = -mass(i) * mass(j) / rmag 
         end if
      end do
      !$omp end parallel do

      pe = pe + sum(pepl(:))

      ! Potential energy from the oblateness term
      if (j2rp2 /= 0.0_DP) then
         !$omp simd 
         do i = 2, npl
            irh(i) = 1.0_DP / norm2(xh(:,i))
         end do
         call obl_pot(npl, swiftest_plA, j2rp2, j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
         pe = pe + oblpot
      end if

      te = ke + pe
   end associate
   return

end subroutine symba_energy_eucl