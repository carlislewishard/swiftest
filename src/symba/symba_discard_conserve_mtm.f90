subroutine symba_discard_conserve_mtm(param, swiftest_plA, ipl, lescape)
   !! author: David A. Minton
   !! 
   !! Conserves system momentum when a body is lost from the system or collides with central body

   use swiftest
   use module_interfaces, except_this_one => symba_discard_conserve_mtm
   implicit none
   ! Arguments
   type(swiftest_pl), intent(inout) :: swiftest_plA
   type(user_input_parameters), intent(inout) :: param
   integer(I4B), intent(in)    :: ipl
   logical, intent(in)         :: lescape
   ! Internals
   real(DP), dimension(NDIM) :: Lpl, Lcb, xcom, vcom
   reaL(DP)                  :: pe, ke

   associate(npl => swiftest_plA%nbody, xb => swiftest_plA%xb, vb => swiftest_plA%vb, &
      rot => swiftest_plA%rot, Ip => swiftest_plA%Ip, radius => swiftest_plA%radius, mass => swiftest_plA%mass, &
      statipl => swiftest_plA%status(ipl), xbpl => swiftest_plA%xb(:,ipl), xbcb => swiftest_plA%xb(:,1))

      xcom(:) = (mass(ipl) * xb(:, ipl) + mass(1) * xb(:,1)) / (mass(1) + mass(ipl))
      vcom(:) = (mass(ipl) * vb(:, ipl) + mass(1) * vb(:,1)) / (mass(1) + mass(ipl))

      call util_crossproduct(xb(:,ipl) - xcom(:), vb(:,ipl) - vcom(:), Lpl)
      Lpl(:) = mass(ipl) * (Lpl(:) + radius(ipl)**2 * Ip(3,ipl) * rot(:, ipl))

      call util_crossproduct(xb(:, 1) - xcom(:), vb(:, 1) - vcom(:), Lcb)
      Lcb(:) = mass(1) * Lcb(:) 

      ! Add planet mass to central body accumulator
      if (lescape) then
         ! Add the potential and kinetic energy of the lost body to the records
         pe = -mass(1) * mass(ipl) / norm2(xb(:, ipl) - xb(:, 1))
         ke = 0.5_DP * mass(ipl) * dot_product(vb(:, ipl), vb(:, ipl))
      else
         ! Add the potential energy of the lost body to the records
         pe = -mass(1) * mass(ipl) / norm2(xb(:, ipl) - xb(:, 1))
         ke = 0.0_DP
         swiftest_plA%dMcb = swiftest_plA%dMcb + mass(ipl)
         swiftest_plA%dRcb = swiftest_plA%dRcb + 1.0_DP / 3.0_DP * (radius(ipl) / radius(1))**3 - 2.0_DP / 9.0_DP * (radius(ipl) / radius(1))**6
         ! Update mass of central body to be consistent with its total mass
         mass(1) = swiftest_plA%Mcb_initial + swiftest_plA%dMcb
         radius(1) = swiftest_plA%Rcb_initial + swiftest_plA%dRcb
         param%rmin = radius(1)
      end if
      swiftest_plA%Ecollisions  = swiftest_plA%Ecollisions - (ke + pe)
   
      ! Add planet angular momentum to central body accumulator
      swiftest_plA%dLcb(:) = Lpl(:) + Lcb(:) + swiftest_plA%dLcb(:)

      ! Update rotation of central body to by consistent with its angular momentum 
      swiftest_plA%rot(:,1) = (swiftest_plA%Lcb_initial(:) + swiftest_plA%dLcb(:)) / (Ip(3, 1) * mass(1) * radius(1)**2)        
      
      ! Update position and velocity of central body
      xb(:, 1) = xcom(:)
      vb(:, 1) = vcom(:)

      ! Update the heliocentric coordinates of everything else
      call coord_b2h(npl, swiftest_plA)
   end associate
   return

end subroutine symba_discard_conserve_mtm