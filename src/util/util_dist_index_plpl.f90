subroutine util_dist_index_plpl(npl, nplm, symba_plA) 
   !! author: Jakob R. Elliott and David A. Minton
   !!
   !! Turns i,j indices into k index for use in the Euclidean distance matrix
   !!  
! modules
   use swiftest
   use swiftest_globals
   use module_symba
   use swiftest_data_structures
   use module_interfaces, EXCEPT_THIS_ONE => util_dist_index_plpl
   implicit none

! arguments
   integer(I4B), intent(in)  :: npl, nplm
   type(symba_pl), intent(inout) :: symba_plA

! internals
   integer(I8B)          :: i,j, counter, npl8, nplm8

   npl8 = int(npl, kind = I8B)
   nplm8 = int(nplm, kind = I8B)
! executable code
   symba_plA%num_plpl_comparisons = ((npl8 - 1_I8B) * (npl8 - 2_I8B) / 2_I8B) - & ! number of entries in a strict lower triangle, nplm x npl, minus first column
               ((npl8 - nplm8 - 1_I8B) * ((npl8 - nplm8 - 1_I8B) + 1_I8B) / 2_I8B)
   if (allocated(symba_plA%k_plpl)) deallocate(symba_plA%k_plpl) 
   allocate(symba_plA%k_plpl(2, symba_plA%num_plpl_comparisons))
   if (allocated(symba_plA%l_plpl_encounter)) deallocate(symba_plA%l_plpl_encounter) 
   allocate(symba_plA%l_plpl_encounter(symba_plA%num_plpl_comparisons))
   ! this is a 'fancier' code, but so far i think it runs slower
   ! so leaving it in, but commenting it out
   ! i think it's because of the 'mod' call, but i haven't profiled it yet
   ! don't forget to uncomment the 'k' declaration up top!
   ! allocate(k(num_comparisons))

   ! m = ceiling(sqrt(2. * num_comparisons))

   ! k = (/(i, i=1,num_comparisons, 1)/)

   ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
   ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2

   ! brute force the index creation

   do i = 2,nplm
      counter = (i - 2) * npl - i * (i - 1) / 2 + 2
      do j = i+1,npl
         symba_plA%k_plpl(1, counter) = i
         symba_plA%k_plpl(2, counter) = j
         counter = counter + 1
      enddo
   enddo

   return

end subroutine util_dist_index_plpl