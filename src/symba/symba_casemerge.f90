subroutine symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list,  & 
   symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Merge planets.
   !! 
   !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
   !!
   !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, except_this_one => symba_casemerge
   implicit none

   integer(I4B), intent(in)                :: index_enc, nplplenc
   integer(I4B), intent(inout)             :: nmergeadd, nmergesub
   real(DP), intent(in)                    :: t
   real(DP), intent(inout)                 :: m1, m2, rad1, rad2
   real(DP), dimension(:), intent(inout)   :: x1, x2, v1, v2
   type(symba_plplenc), intent(inout)      :: plplenc_list
   type(symba_merger), intent(inout)       :: mergeadd_list, mergesub_list
   type(symba_pl), intent(inout)           :: symba_plA
   integer(I4B), dimension(:), intent(in)  :: array_index1_child, array_index2_child

   integer(I4B)                            :: i, j, k, stat1, stat2, index1, index2, indexchild
   integer(I4B)                            :: index1_child, index2_child, index1_parent, index2_parent
   integer(I4B)                            :: name1, name2, nchild1, nchild2
   real(DP)                                :: mtot, Mcb
   real(DP), dimension(NDIM)               :: xnew, vnew, vbs
   integer(I4B), dimension(:), allocatable :: array_keep_child, array_rm_child

   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_plA%index_parent(index1)
   index2_parent = symba_plA%index_parent(index2)
   name1 = symba_plA%helio%swiftest%name(index1)
   name2 = symba_plA%helio%swiftest%name(index2)
   stat1 = symba_plA%helio%swiftest%status(index1)
   stat2 = symba_plA%helio%swiftest%status(index2)
   nchild1 = symba_plA%nchild(index1_parent)
   nchild2 = symba_plA%nchild(index2_parent)
   vbs = symba_plA%helio%swiftest%vb(:, 1)
   Mcb = symba_plA%helio%swiftest%mass(1)

   ! The new position and velocity determined from a perfectly inelastic collision
   mtot = m1 + m2
   xnew(:) = (m1 * x1(:) + m2 * x2(:)) / mtot
   vnew(:) = (m1 * v1(:) + m2 * v2(:)) / mtot
   
   write(*, *) "Merging particles ", name1, " and ", name2, " at time t = ",t

   do k = 1, nplplenc
      if (plplenc_list%status(k) == active) then
         do i = 0, nchild1
            if (i == 0) then 
               index1_child = index1_parent
            else
               index1_child = array_index1_child(i)
            end if 
            do j = 0, nchild2
               if (j == 0) then
                  index2_child = index2_parent
               else
                  index2_child = array_index2_child(j)
               end if

               if ((index1_child == plplenc_list%index1(k)) .and. (index2_child == plplenc_list%index2(k))) then
                  plplenc_list%status(k) = merged
               else if ((index1_child == plplenc_list%index2(k)) .and. (index2_child == plplenc_list%index1(k))) then
                  plplenc_list%status(k) = merged
               end if
            end do
         end do
      end if
   end do

   symba_plA%helio%swiftest%xh(:,index1_parent) = xnew(:)
   symba_plA%helio%swiftest%vb(:,index1_parent) = vnew(:)

   symba_plA%helio%swiftest%xh(:,index2_parent) = xnew(:)
   symba_plA%helio%swiftest%vb(:,index2_parent) = vnew(:) 

   ! the children of parent one are the children we are keeping
   if (nchild1 > 0) then
      allocate(array_keep_child(nchild1))
      array_keep_child(:) = symba_plA%index_child(1:nchild1, index1_parent)
      ! go through the children of the kept parent and add those children to the array of kept children
      do i = 1, nchild1
         indexchild = array_keep_child(i)
         symba_plA%helio%swiftest%xh(:, indexchild) = xnew(:) 
         symba_plA%helio%swiftest%vb(:, indexchild) = vnew(:)
      end do
   end if

   ! the removed parent is assigned as a new child to the list of children of the kept parent
   ! gives kept parent a new child 
   symba_plA%index_child((nchild1 + 1), index1_parent) = index2_parent

   symba_plA%index_parent(index2) = index1_parent
   ! go through the children of the removed parent and add those children to the array of removed children 
   if (nchild2 > 0) then 
      allocate(array_rm_child(nchild2))
      array_rm_child(:) = symba_plA%index_child(1:nchild2, index2_parent)
      ! the parent of the removed parent is assigned to be the kept parent 
      ! gives removed parent a new parent
      do i = 1, nchild2
         symba_plA%index_parent(array_rm_child(i)) = index1_parent
         indexchild = array_rm_child(i)
         symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
         symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
         ! go through the children of the removed parent and add those children to the list of children of the kept parent
         symba_plA%index_child(nchild1 + i + 1, index1_parent) = indexchild
      end do 
   end if
   ! updates the number of children of the kept parent
   symba_plA%nchild(index1_parent) = symba_plA%nchild(index1_parent) + symba_plA%nchild(index2_parent) + 1

   call symba_merger_size_check(mergesub_list, nmergesub + 2)  
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = merged
   mergesub_list%xh(:,nmergesub) = symba_plA%helio%swiftest%xh(:, index1) 
   mergesub_list%vh(:,nmergesub) = symba_plA%helio%swiftest%vh(:, index1) 
   mergesub_list%mass(nmergesub) = m1
   mergesub_list%radius(nmergesub) = rad1
   mergesub_list%nadded(nmergesub) = 1
   mergesub_list%index_ps(nmergesub) = index1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = merged
   mergesub_list%xh(:,nmergesub) = symba_plA%helio%swiftest%xh(:, index2) 
   mergesub_list%vh(:,nmergesub) = symba_plA%helio%swiftest%vh(:, index2)
   mergesub_list%mass(nmergesub) = m2
   mergesub_list%radius(nmergesub) = rad2
   mergesub_list%nadded(nmergesub) = 1
   mergesub_list%index_ps(nmergesub) = index2

   call symba_merger_size_check(mergeadd_list, nmergeadd + 1)  
   nmergeadd = nmergeadd + 1
   if (m2 > m1) then
      mergeadd_list%name(nmergeadd) = name2
      mergeadd_list%status(nmergeadd) = stat2
   else
      mergeadd_list%name(nmergeadd) = name1
      mergeadd_list%status(nmergeadd) = stat1
   end if
   mergeadd_list%ncomp(nmergeadd) = 2
   mergeadd_list%xh(:,nmergeadd) = xnew(:) 
   mergeadd_list%vh(:,nmergeadd) = vnew(:) - vbs(:)
   mergeadd_list%mass(nmergeadd) = mtot
   mergeadd_list%radius(nmergeadd) = (rad1**3 + rad2**3)**(1.0_DP/3.0_DP)

   return 
end subroutine symba_casemerge
