!**********************************************************************************************************************************
!
!  Unit Name   : symba_step
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets and ACTIVE test particles ahead in democratic heliocentric coordinates, descending the recursive
!                branch if necessary to handle possible close encounters
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                t              : time
!                npl            : number of planets
!                ntp            : number of ACTIVE test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of ACTIVE SyMBA test particle structure linked-list
!                param%j2rp2          : J2 * R**2 for the Sun
!                param%j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of ACTIVE SyMBA test particle structure linked-list
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_step(lfirst, param%lextra_force, param%lclose, t, npl, ntp, symba_pl1P, symba_tp1P, param%j2rp2, param%j4rp4,
!                                dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,
!                                mergesub_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_pl.f
!
!**********************************************************************************************************************************
subroutine symba_step(t, dt, param, npl, ntp,symba_plA, symba_tpA,       &
               nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub,&
               mergeadd_list, mergesub_list)
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Step planets and ACTIVE test particles ahead in democratic heliocentric coordinates, descending the recursive
   !! branch if necessary to handle possible close encounters
   !!  
   !! Adapted from Swifter by David E. Kaufmanna symba_step.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_pl.f
! Modules
     use swiftest
     use module_helio
     use module_symba
     use module_swiftestalloc
     use module_interfaces, EXCEPT_THIS_ONE => symba_step
     implicit none

! Arguments
     type(user_input_parameters), intent(inout)       :: param        ! derived type containing user defined parameters 
     integer(I4B), intent(in)                         :: npl, ntp
     integer(I4B), intent(inout)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
     real(DP), intent(in)                             :: t, dt
     type(symba_pl), intent(inout)                    :: symba_plA
     type(symba_tp), intent(inout)                    :: symba_tpA
     type(symba_plplenc), intent(inout)               :: plplenc_list
     type(symba_pltpenc), intent(inout)               :: pltpenc_list
     type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
! Internals
     logical(lgt)              :: lencounter, lvdotr
     integer(I4B)              :: i, j, irec, nplm
     real(DP), dimension(NDIM) :: xr, vr
     logical, save             :: lfirst = .true.

! Executable code
   irec = 0
   call symba_step_reset(npl, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list)
   nplplenc = 0
   npltpenc = 0
   nplm = count(symba_plA%helio%swiftest%mass(1:npl) >= param%mtiny)

   do i = 2, nplm
      !!$omp parallel do schedule(auto) default(private) &
      !!$omp shared(i, npl, nplm, symba_plA, param, dt, irec, plplenc_list, nplplenc)
      do j = i + 1, npl
            xr(:) = symba_plA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
            vr(:) = symba_plA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
            call symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), &
                  symba_plA%helio%swiftest%rhill(j), dt, irec, lencounter, lvdotr)
            if (lencounter) then
               !!$omp critical
               nplplenc = nplplenc + 1
               call symba_plplenc_size_check(plplenc_list, nplplenc)
               plplenc_list%status(nplplenc) = ACTIVE
               plplenc_list%lvdotr(nplplenc) = lvdotr
               plplenc_list%level(nplplenc) = irec
               ! set the first body to be the biggest in the encounter list
               if (symba_plA%helio%swiftest%mass(i) >= symba_plA%helio%swiftest%mass(j)) then
                  plplenc_list%index1(nplplenc) = i
                  plplenc_list%index2(nplplenc) = j
               else
                  plplenc_list%index1(nplplenc) = j
                  plplenc_list%index2(nplplenc) = i
               end if
               symba_plA%nplenc(i) = symba_plA%nplenc(i) + 1
               symba_plA%nplenc(j) = symba_plA%nplenc(j) + 1
               !!$omp end critical
            end if
         end do
      !!$omp end parallel do
         do j = 1, ntp
            xr(:) = symba_tpA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
            vr(:) = symba_tpA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
            call symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), 0.0_DP, dt, irec, lencounter, lvdotr)
            if (lencounter) then
               npltpenc = npltpenc + 1
               call symba_pltpenc_size_check(pltpenc_list, npltpenc)
               symba_plA%ntpenc(i) = symba_plA%ntpenc(i) + 1
               symba_tpA%nplenc(j) = symba_tpA%nplenc(j) + 1
               pltpenc_list%status(npltpenc) = ACTIVE
               pltpenc_list%lvdotr(npltpenc) = lvdotr
               pltpenc_list%level(npltpenc) = irec
               pltpenc_list%indexpl(npltpenc) = i
               pltpenc_list%indextp(npltpenc) = j
            end if
         end do
   end do

     lencounter = ((nplplenc > 0) .or. (npltpenc > 0))
     if (lencounter) then
          call symba_step_interp(t, npl, nplm, ntp, symba_plA, symba_tpA, &
               dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
               nmergesub, mergeadd_list, mergesub_list,  param)
          lfirst = .true.
     else 
          call symba_step_helio(lfirst, param%lextra_force, t, npl, nplm, ntp,&
               symba_plA%helio, symba_tpA%helio, &
               param%j2rp2, param%j4rp4, dt)
     end if

     return

eND SUBROUTINE symba_step