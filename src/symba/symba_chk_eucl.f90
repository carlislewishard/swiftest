subroutine symba_chk_eucl(npl, irec, symba_plA, dt, plplenc_list, nplplenc)
   !! author: Jakob R. Elliott and David A. Minton
   !!
   !! Check for an encounter using the flatten pl-pl data structures
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_chk.f90
   !! Adapted from Hal Levison's Swift routine symba5_chk.f
! modules
   use swiftest
   use swiftest_globals
   use swiftest_data_structures
   use module_helio
   use module_symba
   use module_swiftestalloc
   use module_interfaces, EXCEPT_THIS_ONE => symba_chk_eucl
   !$ use omp_lib
   implicit none

! arguments
   integer(I4B), intent(in)           :: npl, irec
   type(symba_pl), intent(inout)      :: symba_plA
   real(DP), intent(in)               :: dt
   type(symba_plplenc), intent(inout) :: plplenc_list
   integer(I4B), intent(inout)        :: nplplenc

! internals
   logical, dimension(:), allocatable :: loc_lvdotr
   real(DP)   :: r2crit, vdotr, r2, v2, tmin, r2min, term2
   integer(I4B) :: j, nplplenc_new
   integer(I8B) :: k, i
   real(DP), dimension(NDIM):: xr, vr
   integer(I4B), dimension(:,:), allocatable :: ind
   logical, dimension(:), allocatable :: lencounter

! executable code
   associate(xh => symba_plA%helio%swiftest%xh, vh => symba_plA%helio%swiftest%vh, rhill => symba_plA%helio%swiftest%rhill, &
      mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius, &
      num_plpl_comparisons => symba_plA%helio%swiftest%num_plpl_comparisons, k_plpl => symba_plA%helio%swiftest%k_plpl)

      allocate(lencounter(num_plpl_comparisons), loc_lvdotr(num_plpl_comparisons))
      lencounter(:) = .false.

      term2 = rhscale * (rshell**irec)
      nplplenc_new = 0

      do k = 1, num_plpl_comparisons
         associate(ik => k_plpl(1, k), jk => k_plpl(2, k))
            xr(:) = xh(:, jk) - xh(:, ik)
            r2 = dot_product(xr(:), xr(:)) 
            r2crit = ((rhill(jk) + rhill(ik)) * term2)**2
            vr(:) = vh(:, jk) - vh(:, ik)
            vdotr = dot_product(vr(:), xr(:))
            if (r2 < r2crit) then
               lencounter(k) = .true.
               loc_lvdotr(k) = (vdotr < 0.0_DP)
            else
               if (vdotr < 0.0_DP) then
                  v2 = dot_product(vr(:), vr(:))
                  tmin = -vdotr /  v2
                  if (tmin < dt) then
                     r2min = r2 - vdotr * vdotr / v2
                  else
                     r2min = r2 + 2 * vdotr * dt + v2 * dt * dt
                  end if
                  r2min = min(r2min, r2)
                  if (r2min <= r2crit) then
                     lencounter(k) = .true.
                     loc_lvdotr(k) = (vdotr < 0.0_DP)
                  end if
               end if
            end if
         end associate
      end do
      !!$omp end parallel do
      if (any(lencounter(:))) then 
         nplplenc_new = nplplenc + count(lencounter(:))
         call symba_plplenc_size_check(plplenc_list, nplplenc_new)
         plplenc_list%status(nplplenc+1:nplplenc_new) = ACTIVE
         plplenc_list%level(nplplenc+1:nplplenc_new) = irec
         plplenc_list%lvdotr(nplplenc+1:nplplenc_new) = pack(loc_lvdotr(:), lencounter(:))
         plplenc_list%index1(nplplenc+1:nplplenc_new) = pack(k_plpl(1,:), lencounter(:))
         plplenc_list%index2(nplplenc+1:nplplenc_new) = pack(k_plpl(2,:), lencounter(:))
         nplplenc = nplplenc_new
         symba_plA%lencounter(plplenc_list%index1(nplplenc+1:nplplenc_new)) = .true.
         symba_plA%lencounter(plplenc_list%index2(nplplenc+1:nplplenc_new)) = .true.
      end if
      deallocate(loc_lvdotr)
      deallocate(lencounter)
   end associate
   return

   contains

end subroutine symba_chk_eucl