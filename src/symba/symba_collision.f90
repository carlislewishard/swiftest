subroutine symba_collision (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, & 
   Loffset, npl, symba_plA, nplplenc, plplenc_list, mtiny, param)
   !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
   !!
   !! Check for merger between planets in SyMBA. If the user has turned on the FRAGMENTATION feature, it will call the 
   !! symba_regime subroutine to determine what kind of collision will occur.
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine symba_merge_pl.f90
   !!
   !! Adapted from Hal Levison's Swift routine symba5_merge.f
   use swiftest
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => symba_collision
   implicit none

   integer(I4B), intent(in)                   :: index_enc
   integer(I4B), intent(in)                   :: npl, nplplenc
   integer(I4B), intent(inout)                :: nmergeadd, nmergesub
   real(DP), intent(in)                       :: t, dt
   real(DP), intent(inout)                    :: eoffset, Loffset, mtiny
   type(symba_plplenc), intent(inout)         :: plplenc_list
   type(symba_merger), intent(inout)          :: mergeadd_list, mergesub_list
   type(symba_pl), intent(inout)              :: symba_plA
   type(user_input_parameters), intent(inout) :: param

   integer(I4B), parameter                 :: NRES = 3   !! Number of collisional product results
   integer(I4B)                            :: model, i, j, jtarg, jproj
   real(DP), dimension(NRES)               :: mass_res, radius_res, density_res
   real(DP), dimension(NDIM)               :: vbs
   real(DP), dimension(NDIM, NRES)         :: pres, vres
   integer(I4B)                            :: regime, idx_child
   integer(I4B), dimension(2)              :: idx, idx_parent, name, status, nchild 
   real(DP), dimension(2)                  :: radius, mass, density, volume          
   real(DP), dimension(2)                  :: radius_si, mass_si, density_si
   real(DP), dimension(NDIM, 2)            :: x, v, x_si, v_si
   real(DP)                                :: r2, rlim, rlim2, vdotr, tcr2, dt2, a, e, q
   real(DP)                                :: vchild, dentot, Mcb_si
   real(DP)                                :: mmax, mtmp, mtot
   real(DP), dimension(NDIM)               :: xr, vr
   real(DP)                                :: mtiny_si
   logical                                 :: lfrag_add, lmerge
   integer(I4B), dimension(:), allocatable :: array_index1_child, array_index2_child
   real(DP)                                :: mlr, mslr
   integer(I4B)                            :: addi, addf, subi, subf

   ! recalculates vbs 
   call coord_vb2vh(npl, symba_plA%helio%swiftest)
   vbs = symba_plA%helio%swiftest%vb(:, 1)

   lmerge = .false.
   lfrag_add = .false.
   model = 2 ! model 2 is the model for collresolve_resolve (ls12)

   idx(1) = plplenc_list%index1(index_enc)
   idx(2) = plplenc_list%index2(index_enc)

   rlim = symba_plA%helio%swiftest%radius(idx(1)) + symba_plA%helio%swiftest%radius(idx(2))
   xr(:) = symba_plA%helio%swiftest%xh(:, idx(2)) - symba_plA%helio%swiftest%xh(:, idx(1))
   r2 = dot_product(xr(:), xr(:))
   rlim2 = rlim**2
   ! checks if bodies are actively colliding in this time step
   if (r2 <= rlim2) then 
      lfrag_add = .true.
   ! if they are not actively colliding in  this time step, 
   !checks if they are going to collide next time step based on velocities and q
   else 
      vr(:) = symba_plA%helio%swiftest%vb(:, idx(2)) - symba_plA%helio%swiftest%vb(:, idx(1))
      vdotr = dot_product(xr(:), vr(:))
      if (plplenc_list%lvdotr(index_enc) .and. (vdotr > 0.0_DP)) then 
         tcr2 = r2 / dot_product(vr(:), vr(:))
         dt2 = dt**2
         if (tcr2 <= dt2) then
            mtot = symba_plA%helio%swiftest%mass(idx(1)) + symba_plA%helio%swiftest%mass(idx(2))
            call orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
            if (q < rlim) lfrag_add = .true.
         end if
         ! if no collision is going to happen, write as close encounter, not merger
         if (.not. lfrag_add) then
            if (param%encounter_file /= "") then
               name(:)   = symba_plA%helio%swiftest%name(idx(:))
               mass(:)   = symba_plA%helio%swiftest%mass(idx(:))
               radius(:) = symba_plA%helio%swiftest%radius(idx(:))
               do j = 1, 2
                  x(:, j)   = symba_plA%helio%swiftest%xh(:,idx(j)) 
                  v(:, j)   = symba_plA%helio%swiftest%vb(:,idx(j)) - vbs(:)
               end do

               call io_write_encounter(t, name(1), name(2), mass(1), mass(2), radius(1), radius(2), x(:, 1), x(:, 2), &
                  v(:, 1), v(:, 2), param%encounter_file)
            end if
         end if
      end if
   end if
   if (.not. lfrag_add) return ! No collisions, go home.
   
   symba_plA%lmerged(idx(:)) = .true.
   idx_parent(:) = symba_plA%index_parent(idx(:))
   mass(:) = symba_plA%helio%swiftest%mass(idx_parent(:))
   name(:) = symba_plA%helio%swiftest%name(idx_parent(:))
   status(:) = symba_plA%helio%swiftest%status(idx_parent(:))
   radius(:) = symba_plA%helio%swiftest%radius(idx_parent(:))
   volume(:) =  (4.0_DP / 3.0_DP) * pi * radius(:)**3

   nchild(:) = symba_plA%nchild(idx_parent(:)) 
   if (nchild(1) > 0) then
      allocate(array_index1_child, source = symba_plA%index_child(1:nchild(1), idx_parent(1)))
   else 
      allocate(array_index1_child(1))
      array_index1_child(1) = idx_parent(1) 
   end if

   if (nchild(2) > 0) then
      allocate(array_index2_child, source = symba_plA%index_child(1:nchild(2), idx_parent(2)))
   else 
      allocate(array_index2_child(1))
      array_index2_child(1) = idx_parent(2)
   end if

   do j = 1, 2
      x(:, j) = mass(j) * symba_plA%helio%swiftest%xh(:, idx_parent(j))
      v(:, j) = mass(j) * symba_plA%helio%swiftest%vb(:, idx_parent(j))

      mmax = mass(j)
      if (nchild(j) > 0) then
         do i = 1, nchild(j) ! initialize an array of children
            if (j == 1) then
               idx_child = array_index1_child(i)
            else
               idx_child = array_index2_child(i)
            end if
            mtmp = symba_plA%helio%swiftest%mass(idx_child)
            vchild = (4.0_DP / 3.0_DP) * pi * symba_plA%helio%swiftest%radius(idx_child)**3
            volume(j) = volume(j) + vchild
            if (mtmp > mmax) then
               mmax = mtmp
               name(j) = symba_plA%helio%swiftest%name(idx_child)
               status(j) = symba_plA%helio%swiftest%status(idx_child)
            end if
            mass(j) = mass(j) + mtmp
            x(:, j) = x(:, j) + mtmp * symba_plA%helio%swiftest%xh(:, idx_child)
            v(:, j) = v(:, j) + mtmp * symba_plA%helio%swiftest%vb(:, idx_child)
         end do
      end if
      density(j) =  mass(j) / volume(j)
      radius(j) = ((3 * mass(j)) / (density(j) * 4 * pi))**(1.0_DP / 3.0_DP)
      x(:, j) = x(:, j) / mass(j)
      v(:, j) = v(:, j) / mass(j)
   end do
   
   if (param%lfragmentation) then !! If user has enabled this feature, determine the collisional regime and resolve the collision
      ! Convert all quantities to SI units and determine which of the pair is the projectile vs. target before sending them 
      ! to symba_regime
      mass_si(:)    = (mass(:) / GU) * MU2KG 
      radius_si(:)  = radius(:) * DU2M
      x_si(:, :)    = x(:, :) * DU2M
      v_si(:, :)    = v(:, :) * DU2M / TU2S
      density_si(:) = (density(:) / GU) * MU2KG / DU2M**3
      Mcb_si        = symba_plA%helio%swiftest%mass(1) * MU2KG / GU
      mtiny_si      = (mtiny / GU) * MU2KG
   
      mass_res(:) = 0.0_DP
      radius_res(:) = 0.0_DP
      pres(:,:) = 0.0_DP
      vres(:,:) = 0.0_DP

      ! Determine which body is the projectile and which is the target
      if (mass(1) > mass(2)) then
         jtarg = 1
         jproj = 2
      else
         jtarg = 2
         jproj = 1
      end if
      mtot = sum(mass_si(:)) 
      dentot = sum(mass_si(:) * density_si(:)) / mtot 

      call symba_regime(Mcb_si, mass_si(jtarg), mass_si(jproj), radius_si(jtarg), radius_si(jproj), x_si(:, jtarg), x_si(:, jproj),& 
                        v_si(:, jtarg), v_si(:, jproj), density_si(jtarg), density_si(jproj), regime, mlr, mslr, mtiny_si)

      mass_res(1) = min(max(mlr, 0.0_DP), mtot)
      mass_res(2) = min(max(mslr, 0.0_DP), mtot)
      mass_res(3) = min(max(mtot - mlr - mslr, 0.0_DP), mtot)
      density_res(1) = density_si(jtarg)
      density_res(2) = density_si(jproj)
      density_res(3) = dentot

      radius_res(:) = (3 * mass_res(:)  / (4 * pi * density_res(:)))**(1.0_DP / 3.0_DP)
   
      mass_res(:) = (mass_res(:) / MU2KG) * GU
      radius_res(:) = radius_res(:) / DU2M
   else !! When user has *not* enabled FRAGMENTATION, do every collision as a pure merger.
      regime = COLLRESOLVE_REGIME_MERGE
   end if

   if (param%lenergy) then
      ! Save the first index of the newest bodies to be added to the mergeadd/sub lists
      subi = nmergesub + 1 
      addi = nmergeadd + 1
   end if

   call symba_caseresolve(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, vbs, & 
                          symba_plA, nplplenc, plplenc_list, regime, param%plmaxname, param%tpmaxname, &
                          mass_res, radius_res, array_index1_child, array_index2_child, mass(1), mass(2), &
                          radius(1), radius(2), x(:, 1), x(:, 2), v(:, 1), v(:, 2), mtiny, Loffset)


   deallocate(array_index1_child, array_index2_child)
   return

end subroutine symba_collision