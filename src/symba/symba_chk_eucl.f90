subroutine symba_chk_eucl(symba_plA, dt, lvdotr, nplplenc)
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
   use module_interfaces, EXCEPT_THIS_ONE => symba_chk_eucl
   !$ use omp_lib
   implicit none

! arguments
   type(symba_pl), intent(inout)       :: symba_plA
   real(DP), intent(in)                :: dt

   logical(lgt), dimension(:), allocatable, intent(inout) :: lvdotr
   integer(I4B), intent(out)              :: nplplenc

! internals
   logical, dimension(:), allocatable :: loc_lvdotr, ltmp
   ! logical(lgt) :: iflag lvdotr_flag
   real(DP)   :: rcrit, r2crit, vdotr, r2, v2, tmin, r2min, term2, rcritmax, r2critmax
   integer(I4B) :: i, j, npl
   integer(I8B) :: k
   real(DP), dimension(NDIM):: xr, vr

! executable code
   associate(xh => symba_plA%helio%swiftest%xh, vh => symba_plA%helio%swiftest%vh, rhill => symba_plA%helio%swiftest%rhill, &
      mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius)

      nplplenc = 0
      npl = size(symba_plA%helio%swiftest%mass)

      allocate(loc_lvdotr(npl))

      term2 = rhscale*rshell**0

      rcritmax = (symba_plA%helio%swiftest%rhill(2) + symba_plA%helio%swiftest%rhill(3)) * term2
      r2critmax = rcritmax * rcritmax

      do k = 1, symba_plA%num_plpl_comparisons
         associate(ik => symba_plA%k_plpl(1, k), jk => symba_plA%k_plpl(2, k))
            xr(:) = xh(:, jk) - xh(:, ik)
            r2 = dot_product(xr(:), xr(:)) 
            if (r2 < r2critmax) then
               rcrit = (rhill(jk) +rhill(ik)) * term2
               r2crit = rcrit**2 
               vr(:) = vh(:, jk) - vh(:, ik)
               vdotr = dot_product(vr(:), xr(:))
               if (r2 < r2crit) then
                  nplplenc = nplplenc + 1
                  if (nplplenc > size(loc_lvdotr)) then ! Expand array
                     allocate(ltmp(nplplenc * 2))
                     ltmp(1:nplplenc-1) = loc_lvdotr(1:nplplenc-1)
                     call move_alloc(ltmp, loc_lvdotr)
                  end if
                  loc_lvdotr(nplplenc) = (vdotr < 0.0_DP)
                  symba_plA%l_plpl_encounter(k) = .true.
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
                        nplplenc = nplplenc + 1
                        if (nplplenc > size(loc_lvdotr)) then ! Expand array
                           allocate(ltmp(nplplenc * 2))
                           ltmp(1:nplplenc-1) = loc_lvdotr(1:nplplenc-1)
                           call move_alloc(ltmp, loc_lvdotr)
                        end if
                        loc_lvdotr(nplplenc) = (vdotr < 0.0_DP)
                        symba_plA%l_plpl_encounter(k) = .true.
                     end if
                  end if
               end if
            end if
         end associate
      end do
      !!$omp end parallel do
      if (nplplenc > 0) then
         allocate(lvdotr(nplplenc))
         lvdotr(1:nplplenc) = loc_lvdotr(1:nplplenc)
      end if
      deallocate(loc_lvdotr)
   end associate
   return

end subroutine symba_chk_eucl