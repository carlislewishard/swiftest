subroutine symba_discard_conserve_mtm(swiftest_plA, ipl, lescape)
   !! author: David A. Minton
   !! 
   !! Conserves system momentum when a body is lost from the system or collides with central body

   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => symba_discard_conserve_mtm
   implicit none

   integer(I4B), intent(in)    :: ipl
   type(swiftest_pl), intent(inout) :: swiftest_plA
   logical, intent(in)         :: lescape

   real(DP), dimension(NDIM) :: Lpl, Lcb, xcom, vcom

   associate(xb => swiftest_plA%xb, vb => swiftest_plA%vb, rot => swiftest_plA%rot, Ip => swiftest_plA%Ip, rad => swiftest_plA%radius, mass => swiftest_plA%mass, &
      statipl => swiftest_plA%status(ipl), xbpl => swiftest_plA%xb(:,ipl), xbcb => swiftest_plA%xb(:,1))
   
      xcom(:) = (mass(ipl) * xb(:, ipl) + mass(1) * xb(:,1)) / (mass(1) + mass(ipl))
      xcom(:) = (mass(ipl) * vb(:, ipl) + mass(1) * vb(:,1)) / (mass(1) + mass(ipl))

      call util_crossproduct(xb(:,ipl) - xcom(:),vb(:,ipl) - vcom(:), Lpl)
      Lpl(:) = mass(ipl) * (Lpl(:) + rad(ipl)**2 * Ip(3,ipl) * rot(:, ipl))

      call util_crossproduct(xb(:, 1) - xcom(:), vb(:, 1) - vcom(:), Lcb)
      Lcb(:) = mass(1) * Lcb(:) 

      ! Add planet mass to central body accumulator
      if (lescape) then
         swiftest_plA%Mescape = swiftest_plA%Mescape + mass(ipl)
      else
         swiftest_plA%dMcb = swiftest_plA%dMcb + mass(ipl)
         ! Update mass of central body to be consistent with its total mass
         mass(1) = swiftest_plA%Mcb_initial + swiftest_plA%dMcb
      end if
   
      ! Add planet angular momentum to central body accumulator
      swiftest_plA%dLcb(:) = Lpl(:) + Lcb(:) + swiftest_plA%dLcb(:)

      ! Update rotation of central body to by consistent with its angular momentum 
      swiftest_plA%rot(:,1) = (swiftest_plA%Lcb_initial(:) + swiftest_plA%dLcb(:)) / (Ip(3, 1) * mass(1) * rad(1)**2)        
      
      ! Update position and velocity of central body
      !xb(:, 1) = xcom(:)
      !vb(:, 1) = vcom(:)
   end associate
   return

end subroutine symba_discard_conserve_mtm