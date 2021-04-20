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

   real(DP)            :: mass, rad, Ipz, Ipcbz, Mcb, radcb
   real(DP), dimension(NDIM) :: Lpl, Lcb, rot, rotcb, xb, vb, xbcb, vbcb, xcom, vcom

   xb(:) = swiftest_plA%xb(:,ipl)
   vb(:) = swiftest_plA%vb(:,ipl)
   rot(:) = swiftest_plA%rot(:,ipl)
   Ipz = swiftest_plA%Ip(3,ipl)
   rad = swiftest_plA%radius(ipl)
   mass = swiftest_plA%mass(ipl) 

   xbcb(:) = swiftest_plA%xb(:,1)
   vbcb(:) = swiftest_plA%vb(:,1)
   rotcb(:) = swiftest_plA%rot(:,1)
   Ipcbz = swiftest_plA%Ip(3,1)
   radcb = swiftest_plA%radius(1)
   Mcb = swiftest_plA%mass(1)

   xcom(:) = (mass * xb(:) + Mcb * xbcb(:)) / (mass + Mcb)
   vcom(:) = (mass * vb(:) + Mcb * vbcb(:)) / (mass + Mcb)

   call util_crossproduct(xb-xcom,vb-vcom,Lpl)
   Lpl(:) = mass * (Lpl(:) + rad**2 * Ipz * rot(:))

   call util_crossproduct(xbcb-xcom,vbcb-vcom,Lcb)
   Lcb(:) = Mcb * Lcb(:) 

   ! Add planet mass to central body accumulator
   if (lescape) then
      swiftest_plA%Mescape = swiftest_plA%Mescape + mass
   else
      swiftest_plA%dMcb = swiftest_plA%dMcb + mass

      ! Update mass of central body to be consistent with its total mass
      Mcb = swiftest_plA%Mcb_initial + swiftest_plA%dMcb
      swiftest_plA%mass(1) = Mcb
   end if
   
   ! Add planet angular momentum to central body accumulator
   swiftest_plA%dLcb(:) = Lpl(:) + Lcb(:) + swiftest_plA%dLcb(:)

   ! Update rotation of central body to by consistent with its angular momentum 
   swiftest_plA%rot(:,1) = (swiftest_plA%Lcb_initial(:) + swiftest_plA%dLcb(:)) / (Ipcbz * Mcb * radcb**2)        
      
   ! Update position and velocity of central body
   swiftest_plA%xb(:,1) = xcom(:)
   swiftest_plA%vb(:,1) = vcom(:)

   return

end subroutine symba_discard_conserve_mtm