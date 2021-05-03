subroutine symba_merge_pl(t, dt, index_enc, nmergesub, mergesub_list, npl, symba_plA, plplenc_list, param)
   !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
   !!
   !! Check for a collision between planets in SyMBA. This subroutine only flags collising bodies, but does not resolve the outcome.
   !! The collisional outcomes are determined after the step is done in symba_collision
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine symba_merge_pl.f90
   !!
   !! Adapted from Hal Levison's Swift routine symba5_merge.f
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc
   use module_interfaces, EXCEPT_THIS_ONE => symba_merge_pl
   implicit none

   integer(I4B), intent(in)                   :: index_enc
   integer(I4B), intent(in)                   :: npl
   integer(I4B), intent(inout)                :: nmergesub
   real(DP), intent(in)                       :: t, dt
   type(symba_plplenc), intent(inout)         :: plplenc_list
   type(symba_merger), intent(inout)          :: mergesub_list
   type(symba_pl), intent(inout)              :: symba_plA
   type(user_input_parameters), intent(inout) :: param

   integer(I4B)                               :: i, j, k, index_parent, index_child, p1, p2
   integer(I4B)                               :: nchild_inherit, nchild_orig, nchild_new
   real(DP), dimension(NDIM)                  :: vbs
   integer(I4B), dimension(2)                 :: idx, name
   real(DP), dimension(2)                     :: radius, mass 
   real(DP)                                   :: r2, rlim, rlim2, vdotr, tcr2, dt2, a, e, q
   real(DP), dimension(NDIM)                  :: xr, vr
   logical                                    :: lcollision
   integer(I4B), dimension(:), allocatable    :: temp
   real(DP)                                   :: mtot


   ! The plplenc_list is populated such that the most massive of the two bodies is index1
   idx(1) = plplenc_list%index1(index_enc)
   idx(2) = plplenc_list%index2(index_enc)

   rlim = symba_plA%helio%swiftest%radius(idx(1)) + symba_plA%helio%swiftest%radius(idx(2))
   xr(:) = symba_plA%helio%swiftest%xh(:, idx(2)) - symba_plA%helio%swiftest%xh(:, idx(1))
   r2 = dot_product(xr(:), xr(:))
   rlim2 = rlim**2

   lcollision = .false.
   ! checks if bodies are actively colliding in this time step
   if (r2 <= rlim2) then 
      lcollision = .true.
      ! if they are not actively colliding in  this time step, checks if they are going to collide next time step based on velocities and q
   else 
      vr(:) = symba_plA%helio%swiftest%vb(:, idx(2)) - symba_plA%helio%swiftest%vb(:, idx(1))
      vdotr = dot_product(xr(:), vr(:))
      if (plplenc_list%lvdotr(index_enc) .and. (vdotr > 0.0_DP)) then 
         tcr2 = r2 / dot_product(vr(:), vr(:))
         dt2 = dt**2
         if (tcr2 <= dt2) then
            mtot = sum(symba_plA%helio%swiftest%mass(idx(:)))
            call orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
            if (q < rlim) then 
               lcollision = .true.
            end if
         end if

      end if
   end if

   if (lcollision) then
      plplenc_list%status(index_enc) = COLLISION
      plplenc_list%xh1(:,index_enc) = symba_plA%helio%swiftest%xh(:,idx(1)) 
      plplenc_list%vb1(:,index_enc) = symba_plA%helio%swiftest%vb(:,idx(1)) 
      plplenc_list%xh2(:,index_enc) = symba_plA%helio%swiftest%xh(:,idx(2))
      plplenc_list%vb2(:,index_enc) = symba_plA%helio%swiftest%vb(:,idx(2))   
      if (any(symba_plA%lcollision(idx(:)))) then ! One of these bodies has been involved in a collision before, so we need to build a parent/child relationship for them
         p1 = symba_plA%kin(idx(1))%parent
         p2 = symba_plA%kin(idx(2))%parent
         if (p1 == p2) return ! This is a collision between to children of a shared parent. We will ignore it.

         if (symba_plA%helio%swiftest%mass(p1) > symba_plA%helio%swiftest%mass(p2)) then
            index_parent = p1
            index_child = p2
         else
            index_parent = p2
            index_child = p1
         end if

         ! Expand the child array (or create it if necessary) and copy over the previous lists of children
         nchild_orig = symba_plA%kin(index_parent)%nchild
         nchild_inherit = symba_plA%kin(index_child)%nchild
         nchild_new = nchild_orig + nchild_inherit + 1
         allocate(temp(nchild_new))

         if (nchild_orig > 0) temp(1:nchild_orig) = symba_plA%kin(index_parent)%child(1:nchild_orig)
         ! Find out if the child body has any children of its own. The new parent wil inherit these children
         if (nchild_inherit > 0) then
            temp(nchild_orig+1:nchild_orig+nchild_inherit) = symba_plA%kin(index_child)%child(1:nchild_inherit)
            do i = 1, nchild_inherit
               j = symba_plA%kin(index_child)%child(i)
               ! Set the childrens' parent to the new parent
               symba_plA%kin(j)%parent = index_parent
            end do
         end if
         if (allocated(symba_plA%kin(index_child)%child)) deallocate(symba_plA%kin(index_child)%child)
         symba_plA%kin(index_child)%nchild = 0
         ! Add the new child to its parent
         symba_plA%kin(index_child)%parent = index_parent
         temp(nchild_new) = index_child
         ! Save the new child array to the parent
         symba_plA%kin(index_parent)%nchild = nchild_new
         call move_alloc(from=temp, to=symba_plA%kin(index_parent)%child)
      end if

      do i = 1, 2
         if (.not.symba_plA%lcollision(idx(i))) then
            ! This is a new collision for this body, so save it to the subtraction list
            symba_plA%lcollision(idx(i)) = .true.
            call symba_merger_size_check(mergesub_list, nmergesub + 1)  
            nmergesub = nmergesub + 1
            mergesub_list%status(nmergesub) = MERGED
            mergesub_list%id(nmergesub) = symba_plA%helio%swiftest%id(idx(i)) 
            mergesub_list%xb(:,nmergesub) = symba_plA%helio%swiftest%xb(:, idx(i)) 
            mergesub_list%vb(:,nmergesub) = symba_plA%helio%swiftest%vb(:, idx(i)) 
            mergesub_list%mass(nmergesub) = symba_plA%helio%swiftest%mass(idx(i))
            mergesub_list%radius(nmergesub) = symba_plA%helio%swiftest%radius(idx(i))
            mergesub_list%nadded(nmergesub) = 1
            mergesub_list%index_ps(nmergesub) = idx(i)
         end if
      end do

   !else ! Not going to collide, so flag this as a close encounter
      !if (param%encounter_file /= "") then
      !   id(:)   = symba_plA%helio%swiftest%id(idx(:))
      !   mass(:)   = symba_plA%helio%swiftest%mass(idx(:))
      !   radius(:) = symba_plA%helio%swiftest%radius(idx(:))
      !end if
   end if

end subroutine symba_merge_pl