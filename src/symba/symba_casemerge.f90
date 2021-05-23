function symba_casemerge (symba_plA, family, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, param) result(status)
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
   ! Arguments
   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(:), intent(in)    :: family
   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:),   intent(in)      :: mass, radius
   type(user_input_parameters),intent(inout) :: param
   ! Result
   integer(I4B)                              :: status
   ! Internals
   integer(I4B)                              :: i, j, mergeid, ibiggest, nfamily
   real(DP)                                  :: mass_new, radius_new, volume_new, pe
   real(DP), dimension(NDIM)                 :: xcom, vcom, xc, vc, xcrossv
   real(DP), dimension(2)                    :: vol
   real(DP), dimension(NDIM)                 :: L_orb_old, L_spin_old
   real(DP), dimension(NDIM)                 :: L_spin_new, rot_new, Ip_new

   associate(Mpl => symba_plA%helio%swiftest%mass, id => symba_plA%helio%swiftest%id, info => symba_plA%helio%swiftest%info, &
             xb => symba_plA%helio%swiftest%xb, Euntracked => symba_plA%helio%swiftest%Euntracked, Ecollisions => symba_plA%helio%swiftest%Ecollisions)

      status = MERGED

      mass_new = sum(mass(:))

      ! The merged body's name will be that of the largest of the two parents 
      ibiggest = maxloc(Mpl(family(:)), dim=1)
      mergeid = id(family(ibiggest))

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
      do i = 1, 2
         xc(:) = x(:, i) - xcom(:)
         vc(:) = v(:, i) - vcom(:)
         call utiL_crossproduct(xc(:), vc(:), xcrossv(:))
         L_orb_old(:) = L_orb_old(:) + mass(i) * xcrossv(:)
      end do

      ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
      L_spin_new(:) = L_orb_old(:) + L_spin_old(:) 

      ! Assume prinicpal axis rotation on 3rd Ip axis
      rot_new(:) = L_spin_new(:) / (Ip_new(3) * mass_new * radius_new**2)

      ! Keep track of the component of potential energy due to the pre-impact family for book-keeping
      nfamily = size(family(:))
      pe = 0.0_DP
      do j = 1, nfamily
         do i = j + 1, nfamily
            pe = pe - Mpl(i) * Mpl(j) / norm2(xb(:, i) - xb(:, j))
         end do
      end do
      Ecollisions  = Ecollisions + pe 
      Euntracked = Euntracked - pe 

      ! Populate the list of new bodies
      call symba_merger_size_check(mergeadd_list, nmergeadd + 1)  
      nmergeadd = nmergeadd + 1
      mergeadd_list%id(nmergeadd) = mergeid
      mergeadd_list%status(nmergeadd) = status
      mergeadd_list%xb(:,nmergeadd) = xcom(:)
      mergeadd_list%vb(:,nmergeadd) = vcom(:)
      mergeadd_list%mass(nmergeadd) = mass_new
      mergeadd_list%radius(nmergeadd) = radius_new
      mergeadd_list%Ip(:,nmergeadd) = Ip_new(:)
      mergeadd_list%rot(:,nmergeadd) = rot_new(:)
      mergeadd_list%info(nmergeadd) = info(family(ibiggest)) 

   end associate

   return 
end function symba_casemerge
