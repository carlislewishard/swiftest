subroutine tides_getacch(t, npl, nplm, symba_plA)
	!! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
	!!
	!! Calculated tidal torques from central body to any planet and from any planet to central body
	!! planet - planet interactions are considered negligeable
	!! Adapted from Mercury-T code from Bolmont et al. 2015
   use swiftest
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => tides_getacch
   implicit none

   real(DP), intent(in)                   :: t
   integer(I4B), intent(in)               :: npl, nplm
   type(symba_pl)                         :: symba_plA
   type(swiftest_parameters),intent(inout) :: param

   do i = 2, nplm
   	dx(:) = symba_plA%helio%swiftest%xh(:,i) - symba_plA%helio%swiftest%xh(:,1)
   	rj2 = DOT_PRODUCT(dx(:), dx(:))
   	ej = dx/sqrt(rj2)
   	vj = symba_plA%helio%swiftest%vh(:,i) - symba_plA%helio%swiftest%vh(:,1)
   	rotj = symba_plA%helio%swiftest%rot(:,i)
   	rot_central = symba_plA%helio%swiftest%rot(:,1)
      Ftr = 
      Pto = 
      Pto_central =  !Eq 5 Bolmont et al. 2015 
      F_tot(:,i) = (Ftr + (Pto + Pto_central) * dot_product(vj, ej)/sqrt(rj2)) * ej + Pto * cross_product((rotj - theta_j), ej) + Pto_central * cross_product((rot_central - theta_j), ej) !Eq 6 Bolmont et al. 2015
      F_central = F_central + F_tot(:,i)
      
   end do 

   do i = 2, nplm
   	acc =  F_tot(:,i) / symba_plA%helio%swiftest%mass(i) + F_central / symba_plA%helio%swiftest%mass(1)
   	symba_plA%helio%ah(:, i) = symba_plA%helio%ah(:, i) + acc ! corrective acceleration: tidal migration 
   end do 

   return

end subroutine tides_getacch 