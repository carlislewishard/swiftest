!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_timestep
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the maximum accurate seed_timestep for the seed growth and migration Runge-Kutta method
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_seed_timestep(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_seed_timestep(swifter_pl1P,ring,seeds,dtin) result(dtout)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_timestep
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(in)       :: ring
      type(ringmoons_seeds), intent(in)      :: seeds
      real(DP), intent(in)                   :: dtin
      real(DP)                               :: dtout

! Internals
      integer(I4B)                           :: i,nfz,rbin
      real(DP)                               :: dGm_max,da_max,sigavg,sig_max,nu_max,n
      real(DP)                               :: torque_term,e,inc
      real(DP),dimension(seeds%N)            :: mdot
      real(DP),dimension(0:ring%N+1)         :: Tlind
      type(ringmoons_ring)                      :: iring
      type(ringmoons_seeds)                     :: iseeds
      

! Executable code

      e = 0.0_DP
      inc = 0.0_DP
      dtout = dtin

      iring%N = ring%N
      iseeds%N = seeds%N
      call ringmoons_allocate(iring,iseeds)
      iring = ring 
      iseeds = seeds

      dGm_max = -1._DP
      
      do i = 1,iseeds%N
         !nfz = seeds%fz_bin_outer(i) - seeds%fz_bin_inner(i) + 1
         rbin = iseeds%rbin(i)
         mdot(i) = ringmoons_seed_dMdt(iring,swifter_pl1P%mass,iring%Gsigma(rbin),iseeds%Gm(i),iseeds%a(i))
         Tlind(:) = ringmoons_lindblad_torque(swifter_pl1P,iring,iseeds%Gm(i),iseeds%a(i),e,inc)
         n = sqrt((swifter_pl1P%mass + iseeds%Gm(i)) / iseeds%a(i)**3)
         iseeds%Ttide(i) = ringmoons_tidal_torque(swifter_pl1P,iseeds%Gm(i),n,iseeds%a(i),e,inc) 
         iseeds%Torque(i) = iseeds%Ttide(i) - sum(Tlind(:)) 
         iseeds%Torque(i) = iseeds%Torque(i) + mdot(i) * iring%Iz(rbin) * iring%w(rbin)
         
         dGm_max = max(dGm_max,mdot(i) / iseeds%Gm(i))
         
      end do
      if (dGm_max > 0.0_DP) then
         dtout = min(dtout,RK_FACTOR / dGm_max)  ! smallest seed_timestep for the seed growth equation 
      end if

      ! Now aim for seed migration accuracy
         !write(*,*) 'migration'
      da_max = maxval(abs(ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(1:seeds%N),iseeds%a(1:seeds%N),&
                          iseeds%Torque(1:seeds%N),mdot(1:seeds%N))) / iseeds%a(1:seeds%N)) 
                     
      if (da_max > 0.0_DP) then
         dtout = min(dtout, RK_FACTOR / da_max) ! smallest step that keeps the body within a approximately single bin size 
      end if

      call ringmoons_deallocate(iring,iseeds)
      return
end function ringmoons_seed_timestep
