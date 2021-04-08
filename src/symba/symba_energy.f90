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
   use module_interfaces, EXCEPT_THIS_ONE => symba_energy
   implicit none

! arguments
   integer(I4B), intent(in)         :: npl
   real(DP), intent(in)             :: j2rp2, j4rp4
   real(DP), intent(out)            :: ke, pe, te, msys
   real(DP), dimension(:), intent(out) :: Ltot
   type(swiftest_pl), intent(inout)     :: swiftest_plA

! internals
   integer(I4B)              :: i, j
   real(DP)                  :: rmag, v2, rot2, oblpot
   real(DP), dimension(NDIM) :: h, dx
   real(DP), dimension(npl)  :: irh


! executable code

   call coord_h2b(npl, swiftest_plA, msys)

   associate(xb => swiftest_plA%xb, vb => swiftest_plA%vb, mass => swiftest_plA%mass, radius => swiftest_plA%radius, &
             Ip => swiftest_plA%Ip, rot => swiftest_plA%rot, xh => swiftest_plA%xh, status => swiftest_plA%status)
      Ltot = 0.0_DP
      ke = 0.0_DP
      !!$omp simd private(v2,rot2,h) reduction(+:ke,Ltot)
      do i = 1, npl
         if (status(i) /=ACTIVE) cycle
         v2 = dot_product(vb(:,i), vb(:,i))
         rot2 = dot_product(rot(:,i), rot(:,i))
         h(1) = xb(2,i) * vb(3,i) - xb(3,i) * vb(2,i)
         h(2) = xb(3,i) * vb(1,i) - xb(1,i) * vb(3,i)
         h(3) = xb(1,i) * vb(2,i) - xb(2,i) * vb(1,i)

         ! Angular momentum from orbit and spin
         Ltot(:) = Ltot(:) + mass(i) * (h(:) + Ip(3,i) * radius(i)**2 * rot(:,i))

         ! Kinetic energy from orbit and spin
         ke = ke + 0.5_DP * mass(i) * (v2 + Ip(3,i) * radius(i)**2 * rot2)
      end do

      pe = 0.0_DP
      !!$omp parallel do default(private) &
      !!$omp shared (xb, mass, npl) &
      !!$omp reduction (-:pe)
      do i = 1, npl - 1
         if (status(i) /=ACTIVE) cycle
         do j = i + 1, npl
            if (status(j) /=ACTIVE) cycle
            dx(:) = xb(:, j) - xb(:, i) 
            rmag = norm2(dx(:)) 
            if (rmag > tiny(rmag)) pe = pe - mass(i) * mass(j) / rmag 
         end do
      end do
      !!$omp end parallel do

      if (j2rp2 /= 0.0_DP) then
         !!$omp simd private(rmag)
         do i = 2, npl
            if (status(i) /=ACTIVE) cycle
            rmag = norm2(xh(:,i))
            irh(i) = 1.0_DP / rmag
         end do
         call obl_pot(npl, swiftest_plA, j2rp2, j4rp4, xh(:,:), irh, oblpot)
         pe = pe + oblpot
      end if
   end associate

   te = ke + pe

   return

end subroutine symba_energy
