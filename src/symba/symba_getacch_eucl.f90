subroutine symba_getacch_eucl(lextra_force, t, npl, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)
   !! author: Jakob R. Elliott and David A. Minton
   !!
   !! Compute heliocentric accelerations of planets using flattened Euclidean loop.
   !!Accelerations in an encounter are not included here 
   !!  
   !! Adapted from David E. Kaufmann Swifter routine symba_getacch.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f

! modules
   use swiftest
   use swiftest_globals
   use swiftest_data_structures
   use module_helio
   use module_symba
   use module_interfaces, EXCEPT_THIS_ONE => symba_getacch_eucl
   use omp_lib
   implicit none

! arguments
   logical(lgt), intent(in)                    :: lextra_force
   integer(I4B), intent(in)                    :: npl, nplplenc
   real(DP), intent(in)                        :: t, j2rp2, j4rp4
   type(symba_pl), intent(inout)               :: symba_plA
   type(symba_plplenc), intent(inout)          :: plplenc_list


! internals
   integer(I8B)                                 :: i, j, k
   real(DP)                                     :: rji2, irij3, faci, facj, r2, rlim2
   real(DP), dimension(npl)                     :: irh
   real(DP), dimension(:, :), allocatable, save :: aobl
   real(DP)                                     :: dx, dy, dz
!executable code
 
   associate(ah => symba_plA%helio%ah, xh => symba_plA%helio%swiftest%xh, &
              mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius)
      ah(:,:) = 0.0_DP
      do k = 1, symba_plA%num_plpl_comparisons
         if (symba_plA%l_plpl_encounter(k)) cycle
         associate(ik => symba_plA%k_plpl(1, k), jk => symba_plA%k_plpl(2, k))
            dx = xh(1, jk) - xh(1, ik)
            dy = xh(2, jk) - xh(2, ik)
            dz = xh(3, jk) - xh(3, ik)
            rji2 = dx**2 + dy**2 + dz**2
            rlim2 = (radius(ik) + radius(jk))**2
            if (rji2 > rlim2) then
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = mass(ik) * irij3
               facj = mass(jk) * irij3
               ah(1, ik) = ah(1, ik) + facj * dx
               ah(2, ik) = ah(2, ik) + facj * dy
               ah(3, ik) = ah(3, ik) + facj * dz
               ah(1, jk) = ah(1, jk) - faci * dx
               ah(2, jk) = ah(2, jk) - faci * dy
               ah(3, jk) = ah(3, jk) - faci * dz
            end if
         end associate
      end do

      if (j2rp2 /= 0.0_DP) then
         do i = 2, npl
            r2 = dot_product(xh(:,i), xh(:,i))
            irh(i) = 1.0_DP/sqrt(r2)
         end do
         call obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, xh(:,:), irh, aobl)
         do i = 2, npl
            ah(:,i) = ah(:,i) + aobl(:, i) - aobl(:, 1)
         end do
      end if

      if (lextra_force) call symba_user_getacch(t, npl, symba_plA)
   end associate

   return

end subroutine symba_getacch_eucl