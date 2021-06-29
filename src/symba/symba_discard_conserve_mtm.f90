subroutine symba_discard_conserve_mtm(param, swiftest_plA, ipl, lescape_body)
   !! author: David A. Minton
   !! 
   !! Conserves system momentum when a body is lost from the system or collides with central body

   use swiftest
   use module_interfaces, except_this_one => symba_discard_conserve_mtm
   implicit none
   ! Arguments
   type(swiftest_pl), intent(inout) :: swiftest_plA
   type(swiftest_parameters), intent(inout) :: param
   integer(I4B), intent(in)    :: ipl
   logical, intent(in)         :: lescape_body
   ! Internals
   real(DP), dimension(NDIM) :: Lpl, Ltot, Lcb, xcom, vcom
   real(DP)                  :: pe, ke_orbit, ke_spin
   integer(I4B)              :: i, oldstat


   associate(npl => swiftest_plA%nbody, xb => swiftest_plA%xb, vb => swiftest_plA%vb, &
      rot => swiftest_plA%rot, Ip => swiftest_plA%Ip, radius => swiftest_plA%radius, mass => swiftest_plA%mass, &
      Lcb_initial => swiftest_plA%Lcb_initial, dLcb => swiftest_plA%dLcb, Ecollisions => swiftest_plA%Ecollisions, &
      Rcb_initial => swiftest_plA%Rcb_initial, dRcb => swiftest_plA%dRcb, &
      Mcb_initial => swiftest_plA%Mcb_initial, dMcb => swiftest_plA%dMcb, &
      Mescape => swiftest_plA%Mescape, Lescape => swiftest_plA%Lescape, Euntracked => swiftest_plA%Euntracked, &
      status => swiftest_plA%status)

      ! Add the potential and kinetic energy of the lost body to the records
      pe = -mass(1) * mass(ipl) / norm2(xb(:, ipl) - xb(:, 1))
      ke_orbit = 0.5_DP * mass(ipl) * dot_product(vb(:, ipl), vb(:, ipl)) 
      ke_spin  = 0.5_DP * mass(ipl) * radius(ipl)**2 * Ip(3, ipl) * dot_product(rot(:, ipl), rot(:, ipl))

      ! Add the pre-collision ke of the central body to the records
      ! Add planet mass to central body accumulator
      if (lescape_body) then
         Mescape = Mescape + mass(ipl)
         do i = 2, npl
            if (i == ipl) cycle
            pe = pe - mass(i) * mass(ipl) / norm2(xb(:, ipl) - xb(:, i))
         end do

         Ltot(:) = 0.0_DP
         do i = 1, npl
            call util_crossproduct(mass(i) * xb(:,i), vb(:,i), Lpl)
            Ltot(:) = Ltot(:) + Lpl(:)
         end do
         call coord_b2h(npl, swiftest_plA)
         oldstat = status(ipl)
         status(ipl) = INACTIVE
         call coord_h2b(npl, swiftest_plA)
         status(ipl) = oldstat
         do i = 1, npl
            if (i == ipl) cycle
            call util_crossproduct(mass(i) * xb(:,i), vb(:,i), Lpl)
            Ltot(:) = Ltot(:) - Lpl(:) 
         end do 
         Lescape(:) = Lescape(:) + Ltot(:) + mass(ipl) * radius(ipl)**2 * Ip(3, ipl) * rot(:, ipl)

      else
         xcom(:) = (mass(ipl) * xb(:, ipl) + mass(1) * xb(:,1)) / (mass(1) + mass(ipl))
         vcom(:) = (mass(ipl) * vb(:, ipl) + mass(1) * vb(:,1)) / (mass(1) + mass(ipl))
         call util_crossproduct(xb(:,ipl) - xcom(:), vb(:,ipl) - vcom(:), Lpl)
         Lpl(:) = mass(ipl) * (Lpl(:) + radius(ipl)**2 * Ip(3,ipl) * rot(:, ipl))
   
         call util_crossproduct(xb(:, 1) - xcom(:), vb(:, 1) - vcom(:), Lcb)
         Lcb(:) = mass(1) * Lcb(:) 

         ke_orbit = ke_orbit + 0.5_DP * mass(1) * dot_product(vb(:, 1), vb(:, 1)) 
         ke_spin = ke_spin + 0.5_DP * mass(1) * radius(1)**2 * Ip(3, 1) * dot_product(rot(:, 1), rot(:, 1))
         ! Update mass of central body to be consistent with its total mass
         dMcb = dMcb + mass(ipl)
         dRcb = dRcb + 1.0_DP / 3.0_DP * (radius(ipl) / radius(1))**3 - 2.0_DP / 9.0_DP * (radius(ipl) / radius(1))**6
         mass(1) = Mcb_initial + dMcb
         radius(1) = Rcb_initial + dRcb
         param%rmin = radius(1)
         ! Add planet angular momentum to central body accumulator
         dLcb(:) = Lpl(:) + Lcb(:) + dLcb(:)
         ! Update rotation of central body to by consistent with its angular momentum 
         rot(:,1) = (Lcb_initial(:) + dLcb(:)) / (Ip(3, 1) * mass(1) * radius(1)**2)        
         ke_spin  = ke_spin - 0.5_DP * mass(1) * radius(1)**2 * Ip(3, 1) * dot_product(rot(:, 1), rot(:, 1)) 
         xb(:, 1) = xcom(:)
         vb(:, 1) = vcom(:)
         ke_orbit = ke_orbit - 0.5_DP * mass(1) * dot_product(vb(:, 1), vb(:, 1)) 
      end if
      call coord_b2h(npl, swiftest_plA)

      ! We must do this for proper book-keeping, since we can no longer track this body's contribution to energy directly
      if (lescape_body) then
         Ecollisions  = Ecollisions + ke_orbit + ke_spin + pe
         Euntracked  = Euntracked - (ke_orbit + ke_spin + pe)
      else
         Ecollisions  = Ecollisions + pe 
         Euntracked = Euntracked - pe
      end if

   end associate
   return

end subroutine symba_discard_conserve_mtm