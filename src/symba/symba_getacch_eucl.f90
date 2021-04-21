subroutine symba_getacch_eucl(lextra_force, t, npl, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list, &
   num_plpl_comparisons, k_plpl)
   !! author: Jake Elliott and David A. Minton
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
   integer(I8B), intent(in)                    :: num_plpl_comparisons
   real(DP), intent(in)                        :: t, j2rp2, j4rp4
   type(symba_pl), intent(inout)               :: symba_plA
   type(symba_plplenc), intent(inout)          :: plplenc_list
   integer(I4B), dimension(:,:), intent(in)    :: k_plpl


! internals
   integer(I8B)                                 :: i, j, k
   real(DP)                                     :: rji2, irij3, faci, facj, r2, rlim2
   real(DP), dimension(npl)                     :: ahpx, ahpy, ahpz, ahmx, ahmy, ahmz, irh
   real(DP), dimension(nplplenc)                :: aepx, aepy, aepz, aemx, aemy, aemz
   real(DP), dimension(:, :), allocatable, save :: aobl
   real(DP)                                     :: dx, dy, dz
!executable code
 
   ahpx(:) = 0.0_DP
   ahpy(:) = 0.0_DP
   ahpz(:) = 0.0_DP
   ahmx(:) = 0.0_DP
   ahmy(:) = 0.0_DP
   ahmz(:) = 0.0_DP
   
   ! call util_dist_eucl_plpl(symba_plA%helio%swiftest%xh, num_plpl_comparisons, k_plpl, dist_plpl_array) ! does not care about mtiny

! there is floating point arithmetic round off error in this loop
! for now, we will keep it in the serial operation, so we can easily compare
! it to the older swifter versions
   associate(ah => symba_plA%helio%ah, xh => symba_plA%helio%swiftest%xh, &
              mass => symba_plA%helio%swiftest%mass, radius => symba_plA%helio%swiftest%radius)
      !$omp parallel do default(private) schedule(auto) &
      !$omp firstprivate(num_plpl_comparisons, k_plpl, xh, mass, radius) &
      !$omp reduction(+:ahpx, ahpy, ahpz) &
      !$omp reduction(-:ahmx, ahmy, ahmz)
      do k = 1, num_plpl_comparisons
         associate(i => k_plpl(1, k), j=> k_plpl(2, k))
            dx = xh(1,j) - xh(1,i)
            dy = xh(2,j) - xh(2,i)
            dz = xh(3,j) - xh(3,i)
            rlim2 = (radius(i) + radius(j))**2
            rji2 = dx**2 + dy**2 + dz**2
            if (rji2 > rlim2) then !if false, we likely have recent fragments with coincident positions. 
               !  so ignore in this step and let the merge code deal with it in the nex
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = mass(i) * irij3
               facj = mass(j) * irij3
               ahpx(i) = ahpx(i) + facj * dx
               ahpy(i) = ahpy(i) + facj * dy
               ahpz(i) = ahpz(i) + facj * dz
               ahmx(j) = ahmx(j) - faci * dx
               ahmy(j) = ahmy(j) - faci * dy
               ahmz(j) = ahmz(j) - faci * dz
            end if
         end associate
      end do
      !$omp end parallel do
      ah(1,1:npl) = ahpx(:) + ahmx(:)
      ah(2,1:npl) = ahpy(:) + ahmy(:)
      ah(3,1:npl) = ahpz(:) + ahmz(:)

      aepx(:) = 0.0_DP
      aepy(:) = 0.0_DP
      aepz(:) = 0.0_DP
      aemx(:) = 0.0_DP
      aemy(:) = 0.0_DP
      aemz(:) = 0.0_DP

      do i = 1, nplplenc
         associate(index_i => plplenc_list%index1(i), index_j => plplenc_list%index2(i))
            dx = xh(1,index_j) - xh(1,index_i)
            dy = xh(2,index_j) - xh(2,index_i)
            dz = xh(3,index_j) - xh(3,index_i)
            rlim2 = (radius(index_i) + radius(index_j))**2
            rji2 = dx**2 + dy**2 + dz**2
            if (rji2 > rlim2) then !if false, we likely have recent fragments with coincident positions. 
                           !  so ignore in this step and let the merge code deal with it in the next
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = mass(index_i) * irij3
               facj = mass(index_j) * irij3
               ah(1, index_i) = ah(1, index_i) - facj * dx
               ah(2, index_i) = ah(2, index_i) - facj * dy
               ah(3, index_i) = ah(3, index_i) - facj * dz
               ah(1, index_j) = ah(1, index_j) + faci * dx
               ah(2, index_j) = ah(2, index_j) + faci * dy
               ah(3, index_j) = ah(3, index_j) + faci * dz
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