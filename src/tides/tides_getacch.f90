subroutine tides_getacch(t, npl, symba_plA)
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
   integer(I4B), intent(in)               :: npl
   type(symba_pl)                         :: symba_plA
   type(swiftest_parameters),intent(inout) :: param




    return

end subroutine tides_getacch 