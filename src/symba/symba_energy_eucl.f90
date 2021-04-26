subroutine symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke, pe, te, Ltot, msys)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector and kinetic, potential and total system energy
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
   use swiftest
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_energy_eucl
   implicit none

! arguments
   integer(I4B), intent(in)                 :: npl
   real(DP), intent(in)                     :: j2rp2, j4rp4
   real(DP), intent(out)                    :: ke, pe, te, msys
   real(DP), dimension(:), intent(out)      :: Ltot
   type(symba_pl), intent(inout)         :: symba_plA

! internals
   integer(I4B)              :: i, j
   integer(I8B)              :: k
   real(DP)                  :: rmag, v2, rot2, oblpot, hx, hy, hz
   real(DP), dimension(npl)  :: irh, kepl, pecb
   real(DP), dimension(npl) :: Lplx, Lply, Lplz
   real(DP), dimension(symba_plA%helio%swiftest%num_plpl_comparisons) :: pepl 
   logical, dimension(symba_plA%helio%swiftest%num_plpl_comparisons) :: lstatpl
   logical, dimension(npl) :: lstatus

! executable code

   call coord_h2b(npl, symba_plA%helio%swiftest, msys)
   Ltot = 0.0_DP
   ke = 0.0_DP
   associate(xb => symba_plA%helio%swiftest%xb, vb => symba_plA%helio%swiftest%vb, mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius, &
      Ip => symba_plA%helio%swiftest%Ip, rot => symba_plA%helio%swiftest%rot, xh => symba_plA%helio%swiftest%xh, status => symba_plA%helio%swiftest%status)
      kepl(:) = 0.0_DP
      Lplx(:) = 0.0_DP
      Lply(:) = 0.0_DP
      Lplz(:) = 0.0_DP
      lstatus(1:npl) = status(1:npl) /= INACTIVE
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

      ke = 0.5_DP * sum(kepl(1:npl), lstatus(:))
      Ltot(1) = sum(Lplx(1:npl), lstatus(1:npl)) 
      Ltot(2) = sum(Lply(1:npl), lstatus(1:npl)) 
      Ltot(3) = sum(Lplz(1:npl), lstatus(1:npl)) 

      ! Do the central body potential energy component first
      !$omp simd 
      do i = 2, npl
         associate(px => xh(1,i), py => xh(2,i), pz => xh(3,i))
            pecb(i) = - mass(1) * mass(i) / sqrt(px**2 + py**2 + pz**2)
         end associate
      end do

      ! Do the potential energy between pairs of massive bodies
      do k = 1, symba_plA%helio%swiftest%num_plpl_comparisons
         associate(ik => symba_plA%helio%swiftest%k_plpl(1, k), jk => symba_plA%helio%swiftest%k_plpl(2, k))
            pepl(k) = -mass(ik) * mass(jk) / norm2(xb(:, jk) - xb(:, ik)) 
            lstatpl(k) = (lstatus(ik) .and. lstatus(jk))
         end associate
      end do

      pe = sum(pecb(2:npl), lstatus(2:npl)) + sum(pepl(:), lstatpl(:))

      ! Potential energy from the oblateness term
      if (j2rp2 /= 0.0_DP) then
         !$omp simd 
         do i = 2, npl
            irh(i) = 1.0_DP / norm2(xh(:,i))
         end do
         call obl_pot(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, oblpot)
         pe = pe + oblpot
      end if

      te = ke + pe
   end associate
   return

end subroutine symba_energy_eucl