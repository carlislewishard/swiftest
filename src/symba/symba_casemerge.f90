subroutine symba_casemerge (symba_plA, idx_parents, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, xbs, vbs, param)
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
   use module_interfaces, EXCEPT_THIS_ONE => symba_casemerge
   implicit none

   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(2), intent(inout) :: idx_parents
   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   real(DP), dimension(:),   intent(in)      :: mass, radius, xbs, vbs
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   type(user_input_parameters),intent(inout) :: param

   integer(I4B)                            :: j, mergename
   real(DP)                                :: mass_new, radius_new, volume_new
   real(DP), dimension(NDIM)               :: xcom, vcom, xc, vc, xcrossv
   real(DP), dimension(2)                  :: vol
   real(DP), dimension(NDIM)               :: L_orb_old, L_spin_old
   real(DP), dimension(NDIM)               :: L_spin_new, rot_new, Ip_new
   integer(I4B), dimension(:), allocatable :: family
   integer(I4B)                            :: fam_size, istart

   ! Make the list of family members (bodies involved in the collision)
   associate(nchild1 => symba_plA%kin(idx_parents(1))%nchild, nchild2 => symba_plA%kin(idx_parents(2))%nchild)
      fam_size = 2 + nchild1 + nchild2
      allocate(family(fam_size))
      family(1) = idx_parents(1)
      family(2) = idx_parents(2)
      istart = 2 + nchild1

      if (nchild1 > 0) family(3:istart) = symba_plA%kin(idx_parents(1))%child(1:nchild1)
      if (nchild2 > 0) family(istart+1:istart+1+nchild2) = symba_plA%kin(idx_parents(2))%child(1:nchild2)
      mass_new = sum(symba_plA%helio%swiftest%mass(family(:)))
      mergename = symba_plA%helio%swiftest%name(idx_parents(1))
   end associate
   ! Merged body is created at the barycenter of the original bodies
   xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mass_new
   vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mass_new

   ! Get mass weighted mean of Ip and 
   Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mass_new
   vol(:) = 4._DP / 3._DP * PI * radius(:)**3
   volume_new = sum(vol(:))
   radius_new = (3 * volume_new / (4 * PI))**(1._DP / 3._DP)
   
   L_spin_old(:) = L_spin(:,1) + L_spin(:,2)
   L_orb_old(:) = 0.0_DP
   ! Compute orbital angular momentum of pre-impact system
   do j = 1, 2
      xc(:) = x(:, j) - xcom(:)
      vc(:) = v(:, j) - vcom(:)
      call utiL_crossproduct(xc(:), vc(:), xcrossv(:))
      L_orb_old(:) = L_orb_old(:) + mass(j) * xcrossv(:)
   end do

   ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
   L_spin_new(:) = L_orb_old(:) + L_spin_old(:) 

   ! Assume prinicpal axis rotation on 3rd Ip axis
   rot_new(:) = L_spin_new(:) / (Ip_new(3) * mass_new * radius_new**2)

   ! Populate the list of new bodies
   call symba_merger_size_check(mergeadd_list, nmergeadd + 1)  
   nmergeadd = nmergeadd + 1
   mergeadd_list%name(nmergeadd) = mergename
   mergeadd_list%status(nmergeadd) = MERGED
   mergeadd_list%xb(:,nmergeadd) = xcom(:)
   mergeadd_list%vb(:,nmergeadd) = vcom(:)
   mergeadd_list%mass(nmergeadd) = mass_new
   mergeadd_list%radius(nmergeadd) = radius_new
   mergeadd_list%Ip(:,nmergeadd) = Ip_new(:)
   mergeadd_list%rot(:,nmergeadd) = rot_new(:)

   return 
end subroutine symba_casemerge
