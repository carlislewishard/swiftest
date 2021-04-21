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
   real(DP)                  :: rmag, v2, rot2, oblpot, hx, hy, hz
   real(DP), dimension(npl)  :: irh, kepl
   real(DP), dimension(npl) :: Lplx, Lply, Lplz
   real(DP), dimension(num_plpl_comparisons) :: pepl 
   logical, dimension(num_plpl_comparisons) :: lstatpl

! executable code

   call coord_h2b(npl, swiftest_plA, msys)
   Ltot = 0.0_DP
   ke = 0.0_DP
   associate(xb => swiftest_plA%xb, vb => swiftest_plA%vb, mass => swiftest_plA%mass, radius => swiftest_plA%radius, &
      Ip => swiftest_plA%Ip, rot => swiftest_plA%rot, xh => swiftest_plA%xh, status => swiftest_plA%status)
      kepl(:) = 0.0_DP
      Lplx(:) = 0.0_DP
      Lply(:) = 0.0_DP
      Lplz(:) = 0.0_DP
      !$omp simd private(v2, rot2, hx, hy, hz)
      do i = 1, npl
         v2 = dot_product(vb(:,i), vb(:,i))
         rot2 = dot_product(rot(:,i), rot(:,i))
         hx = xb(2,i) * vb(3,i) - xb(3,i) * vb(2,i)
         hy = xb(3,i) * vb(1,i) - xb(1,i) * vb(3,i)
         hz = xb(1,i) * vb(2,i) - xb(2,i) * vb(1,i)

         ! Angular momentum from orbit and spin
         Lplx(i) = mass(i) * (Ip(3,i) * radius(i)**2 * rot(1,i) + hx)
         Lply(i) = mass(i) * (Ip(3,i) * radius(i)**2 * rot(2,i) + hy)
         Lplz(i) = mass(i) * (Ip(3,i) * radius(i)**2 * rot(3,i) + hz)

         ! Kinetic energy from orbit and spin
         kepl(i) = mass(i) * (Ip(3,i) * radius(i)**2 * rot2 + v2)
      end do

      ke = 0.5_DP * sum(kepl(1:npl), status(1:npl) == ACTIVE)
      Ltot(1) = sum(Lplx(1:npl), status(1:npl) == ACTIVE) 
      Ltot(2) = sum(Lply(1:npl), status(1:npl) == ACTIVE) 
      Ltot(3) = sum(Lplz(1:npl), status(1:npl) == ACTIVE) 

      ! Do the central body potential energy component first
      pe = 0.0_DP
      !$omp simd reduction(-:pe)
      do i = 2, npl
         if (status(i) == ACTIVE) pe = pe - mass(1) * mass(i) / norm2(xh(:, i))
      end do

      ! Do the potential energy between pairs of massive bodies
      pepl(:) = 0.0_DP
      !!$omp parallel do default(private) schedule(auto)&
      !!$omp shared(pepl, lstatpl) &
      !!$omp firstprivate(k_plpl, mass, xb, status)
      do k = 1, num_plpl_comparisons
         associate(i => k_plpl(1, k), j=> k_plpl(2, k))
            pepl(k) = -mass(i) * mass(j) / norm2(xb(:, j) - xb(:, i)) 
            lstatpl(k) = (status(i) == ACTIVE) .and. (status(j) == ACTIVE)
         end associate
      end do
      !!$omp end parallel do

      pe = pe + sum(pepl(:), lstatpl(:))

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