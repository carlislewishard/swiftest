subroutine symba_energy(npl, symba_plA, j2rp2, j4rp4, ke_orbit, ke_spin, pe, te, Ltot)
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector, kinetic, potential and total 
   !! system energy, and total system mass
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_energy.f90
   !!  
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
   use swiftest
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_energy
   implicit none

! arguments
   integer(I4B), intent(in)         :: npl
   real(DP), intent(in)             :: j2rp2, j4rp4
   real(DP), intent(out)            :: ke_orbit, ke_spin, pe, te
   real(DP), dimension(:), intent(out) :: Ltot
   type(symba_pl), intent(inout)     :: symba_plA

! internals
   integer(I4B)              :: i, j
   real(DP)                  :: rmag, v2, rot2, oblpot
   real(DP), dimension(NDIM) :: h, dx
   real(DP), dimension(npl)  :: irh, kepl, kespinpl
   real(DP), dimension(NDIM, npl) :: Lpl 
   real(DP), dimension(npl,npl) :: pepl

! executable code

   associate(xb => symba_plA%helio%swiftest%xb, vb => symba_plA%helio%swiftest%vb, mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius, &
             Ip => symba_plA%helio%swiftest%Ip, rot => symba_plA%helio%swiftest%rot, xh => symba_plA%helio%swiftest%xh, status => symba_plA%helio%swiftest%status)
      Lpl(:,:) = 0.0_DP
      kepl(:) = 0.0_DP
      do i = 1, npl
         if (status(i) /= ACTIVE) cycle
         v2 = dot_product(vb(:,i), vb(:,i))
         rot2 = dot_product(rot(:,i), rot(:,i))
         call util_crossproduct(xb(:,i), vb(:,i), h)

         ! Angular momentum from orbit and spin
         Lpl(:, i) = mass(i) * (Ip(3,i) * radius(i)**2 * rot(:,i) + h(:))

         ! Kinetic energy from orbit and spin
         kepl(i) = mass(i) * v2
         kespinpl(i) = mass(i) * Ip(3,i) * radius(i)**2 * rot2 
      end do
      ke_orbit = 0.5_DP * sum(kepl(:))
      ke_spin  = 0.5_DP * sum(kespinpl(:))

      pepl(:,:) = 0.0_DP
      do i = 1, npl - 1
         if (status(i) /= ACTIVE) cycle
         do j = i + 1, npl
            if (status(j) /= ACTIVE) cycle
            dx(:) = xb(:, j) - xb(:, i) 
            rmag = norm2(dx(:)) 
            if (rmag > tiny(rmag)) pepl(i,j) = - mass(i) * mass(j) / rmag 
         end do
      end do
      pe = sum(pack(pepl(:,:), abs(pepl(:,:)) > tiny(pe)))
      Ltot(:) =  sum(Lpl(:,:), dim=2)

      if (j2rp2 /= 0.0_DP) then
         do i = 2, npl
            if (status(i) /= ACTIVE) cycle
            rmag = norm2(xh(:,i))
            irh(i) = 1.0_DP / rmag
         end do
         call obl_pot(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, xh(:,:), irh, oblpot)
         pe = pe + oblpot
      end if
   end associate

   te = ke_orbit + ke_spin + pe

   return

end subroutine symba_energy