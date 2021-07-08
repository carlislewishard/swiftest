submodule (symba) s_symba_helio_getacch
contains
   module procedure symba_helio_getacch
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_helio_getacch.f90
   !! Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
use swiftest
implicit none
   logical , save                 :: lmalloc = .true.
   integer(I4B)                     :: i
   real(DP)                       :: r2
   real(DP), dimension(:), allocatable, save    :: irh
   real(DP), dimension(:, :), allocatable, save :: xh, aobl


! executable code
   if (lflag) then
      do i = 2, npl
         helio_plA%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      end do
      call symba_helio_getacch_int(npl, nplm, helio_plA) 
   end if
   if (param%loblatecb) then
      if (lmalloc) then
         allocate(xh(NDIM, param%nplmax), aobl(NDIM, param%nplmax), irh(param%nplmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         xh(:, i) = helio_plA%xh(:,i)
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(helio_plA, param%j2rp2, param%j4rp4, xh, irh, aobl) 
      do i = 2, npl
         helio_plA%ah(:,i) = helio_plA%ah(:,i) + aobl(:, i) - aobl(:, 1)
      end do
   else
      do i = 2, npl
         helio_plA%ah(:,i) = helio_plA%ah(:,i)
      end do
   end if
   if (lextra_force) call helio_user_getacch(t, npl, helio_plA) 

   return

   end procedure symba_helio_getacch
end submodule s_symba_helio_getacch
