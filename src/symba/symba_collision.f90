subroutine symba_collision (t, symba_plA, nplplenc, plplenc_list, ldiscard, mergeadd_list, nmergeadd, param)
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
   use module_interfaces, EXCEPT_THIS_ONE => symba_collision
   implicit none

   real(DP), intent(in)                      :: t
   integer(I4B), intent(inout)               :: nplplenc, nmergeadd
   type(symba_pl)                            :: symba_pla
   type(symba_plplenc), intent(inout)        :: plplenc_list
   type(symba_merger), intent(inout)         :: mergeadd_list
   logical, intent(inout)                    :: ldiscard
   type(user_input_parameters),intent(inout) :: param

   integer(I4B), parameter                 :: NRES = 3   !! Number of collisional product results
   integer(I4B)                            :: i, j, k, index_enc, jtarg, jproj
   real(DP), dimension(NRES)               :: mass_res
   real(DP), dimension(NDIM)               :: x1_si, v1_si, x2_si, v2_si
   integer(I4B)                            :: regime, idx_child, status
   integer(I4B), dimension(2)              :: idx, idx_parent, nchild, id
   real(DP), dimension(2)                  :: radius, mass, density, volume
   real(DP), dimension(2)                  :: radius_si, mass_si, density_si
   real(DP), dimension(NDIM, 2)            :: x, v, L_spin, Ip
   real(DP)                                :: volchild, dentot, Mcb_si
   real(DP)                                :: mchild, mtot
   real(DP), dimension(NDIM)               :: xc, vc, xcom, vcom, xchild, vchild, xcrossv
   real(DP)                                :: mtiny_si
   real(DP)                                :: mlr, mslr, msys, msys_new, Qloss
   integer(I4B), dimension(:), allocatable :: family
   integer(I4B)                            :: fam_size, istart
   type family_array
      integer(I4B), dimension(:), allocatable :: id
      integer(I4B), dimension(:), allocatable :: idx
   end type family_array
   type(family_array), dimension(2)        :: parent_child_index_array

   ! TESTING
   logical, save   :: lfirst = .true.
   real(DP), save  :: Minitial
   real(DP)        :: Msystem, Madd, Mdiscard


   ! First determine the collisional regime for each colliding pair
   associate(npl => symba_plA%helio%swiftest%nbody, xbpl => symba_plA%helio%swiftest%xb, statpl => symba_plA%helio%swiftest%status)
      if (lfirst) then
         Minitial = sum(symba_plA%helio%swiftest%mass(1:npl))
         lfirst = .false.
      end if 
      ldiscard = any(plplenc_list%status(1:nplplenc) == COLLISION)
      if (.not.ldiscard) return

      ! Recompute central body barycentric velocity
      call coord_h2b(npl, symba_plA%helio%swiftest, msys)

      ! Loop through the list of pl-pl encounters and pick out the collisions
      do index_enc = 1, nplplenc
         if (plplenc_list%status(index_enc) /= COLLISION) cycle ! Not the primary collision for this pl-pl encounter, so skip it
         !if (t > 1.01E+05) then
         !   write(*,*) "We've arrived at a problem"
         !end if

         ! Index values of the original particle pair 
         idx(1) = plplenc_list%index1(index_enc)
         idx(2) = plplenc_list%index2(index_enc)

         if (any(statpl(idx(:)) /= ACTIVE)) cycle ! One of these two bodies is already gone

         ! Index values for the parents of this particle pair
         idx_parent(:) = symba_plA%kin(idx(:))%parent

         nchild(:) = symba_plA%kin(idx_parent(:))%nchild 
         ! If all of these bodies share a parent, but this is still a unique collision, move the last child
         ! out of the parent's position and make it the secondary body
         if (idx_parent(1) == idx_parent(2)) then
            write(*,*) idx_parent(1), "is having a collision with itself for some reason."
            idx_parent(2) = symba_plA%kin(idx_parent(1))%child(nchild(1))
            nchild(1) = nchild(1) - 1
            nchild(2) = 0
            symba_plA%kin(idx_parent(1))%nchild = nchild(1)
         end if

         mass(:) = symba_plA%helio%swiftest%mass(idx_parent(:))
         id(:) = symba_plA%helio%swiftest%id(idx_parent(:))
         radius(:) = symba_plA%helio%swiftest%radius(idx_parent(:))
         volume(:) =  (4.0_DP / 3.0_DP) * PI * radius(:)**3
    
         ! Group together the ids and indexes of each collisional parent and its children
         do j = 1, 2
            allocate(parent_child_index_array(j)%idx(nchild(j)+ 1))
            allocate(parent_child_index_array(j)%id(nchild(j)+ 1))
            associate(idx_arr => parent_child_index_array(j)%idx, &
                      id_arr => parent_child_index_array(j)%id, &
                      ncj => nchild(j), &
                      pl => symba_plA%helio%swiftest, &
                      plkinj => symba_plA%kin(idx_parent(j)))
               idx_arr(1) = idx_parent(j)
               if (ncj > 0) idx_arr(2:ncj + 1) = plkinj%child(1:ncj)
               id_arr(:) = pl%id(idx_arr(:))
            end associate
         end do

         ! Consolidate the groups of collsional parents with any children they may have into a single "family" index array
         fam_size = 2 + sum(nchild(:))
         allocate(family(fam_size))
         family = [parent_child_index_array(1)%idx(:),parent_child_index_array(2)%idx(:)]

         ! Prepare to resolve collisions by setting the status flag for all family members to COLLISION. This will get updated after
         ! we have determined what kind of collision this group will produce.
         where (statpl(family(:)) == ACTIVE) statpl(family(:)) = COLLISION

         ! Find the barycenter of each body along with its children, if it has any
         do j = 1, 2
            x(:, j)  = symba_plA%helio%swiftest%xb(:, idx_parent(j))
            v(:, j)  = symba_plA%helio%swiftest%vb(:, idx_parent(j))
            Ip(:, j) = mass(j) * symba_plA%helio%swiftest%Ip(:, idx_parent(j))
            ! Assume principal axis rotation about axis corresponding to highest moment of inertia (3rd Ip)
            L_spin(:, j)  = Ip(3, j) * radius(j)**2 * symba_plA%helio%swiftest%rot(:, idx_parent(j))
            if (nchild(j) > 0) then
               do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
                  idx_child = parent_child_index_array(j)%idx(i + 1)
                  if (statpl(idx_child) /= COLLISION) cycle
                  mchild = symba_plA%helio%swiftest%mass(idx_child)
                  xchild(:) = symba_plA%helio%swiftest%xb(:, idx_child)
                  vchild(:) = symba_plA%helio%swiftest%vb(:, idx_child)
                  volchild = (4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(idx_child)**3
                  volume(j) = volume(j) + volchild
                  ! Get angular momentum of the child-parent pair and add that to the spin
                  xcom(:) = (mass(j) * x(:,j) + mchild * xchild(:)) / (mass(j) + mchild)
                  vcom(:) = (mass(j) * v(:,j) + mchild * vchild(:)) / (mass(j) + mchild)
                  xc(:) = x(:, j) - xcom(:)
                  vc(:) = v(:, j) - vcom(:)
                  call util_crossproduct(xc(:), vc(:), xcrossv(:))
                  L_spin(:, j) = L_spin(:, j) + mass(j) * xcrossv(:)
                  xc(:) = xchild(:) - xcom(:)
                  vc(:) = vchild(:) - vcom(:)
                  call util_crossproduct(xc(:), vc(:), xcrossv(:))
                  L_spin(:, j) = L_spin(:, j) + mchild * xcrossv(:)

                  ! Add the child's spin
                  L_spin(:, j) = L_spin(:, j) + mchild * symba_plA%helio%swiftest%Ip(3, idx_child) * symba_plA%helio%swiftest%radius(idx_child)**2 * &
                                             symba_plA%helio%swiftest%rot(:, idx_child)

                  ! Merge the child and parent
                  mass(j) = mass(j) + mchild
                  x(:, j) = xcom(:)
                  v(:, j) = vcom(:)
                  Ip(:, j) = Ip(:, j) + mchild * symba_plA%helio%swiftest%Ip(:, idx_child)
               end do
            end if
            density(j) =  mass(j) / volume(j)
            radius(j) = ((3 * mass(j)) / (density(j) * 4 * pi))**(1.0_DP / 3.0_DP)
            Ip(:, j) = Ip(:, j) / mass(j)
         end do

         if (param%lfragmentation) then !! If user has enabled this feature, determine the collisional regime and resolve the collision
            ! Convert all quantities to SI units and determine which of the pair is the projectile vs. target before sending them 
            ! to symba_regime
            if (mass(1) > mass(2)) then
               jtarg = 1
               jproj = 2
            else
               jtarg = 2
               jproj = 1
            end if
            mass_si(:)    = (mass(:) / GU) * MU2KG                            !! The collective mass of the parent and its children
            radius_si(:)  = radius(:) * DU2M                                  !! The collective radius of the parent and its children
            x1_si(:)      = plplenc_list%xh1(:,index_enc) * DU2M              !! The position of the parent from inside the step (at collision)
            v1_si(:)      = plplenc_list%vb1(:,index_enc) * DU2M / TU2S       !! The velocity of the parent from inside the step (at collision)
            x2_si(:)      = plplenc_list%xh2(:,index_enc) * DU2M              !! The position of the parent from inside the step (at collision)
            v2_si(:)      = plplenc_list%vb2(:,index_enc) * DU2M / TU2S       !! The velocity of the parent from inside the step (at collision)
            density_si(:) = (density(:) / GU) * MU2KG / DU2M**3               !! The collective density of the parent and its children
            Mcb_si        = symba_plA%helio%swiftest%mass(1) * MU2KG / GU
            mtiny_si      = (param%mtiny / GU) * MU2KG
         
            mass_res(:) = 0.0_DP
      
            mtot = sum(mass_si(:)) 
            dentot = sum(mass_si(:) * density_si(:)) / mtot 

            !! Use the positions and velocities of the parents from indside the step (at collision) to calculate the collisional regime
            call symba_regime(Mcb_si, mass_si(jtarg), mass_si(jproj), radius_si(jtarg), radius_si(jproj), x1_si(:), x2_si(:),& 
                  v1_si(:), v2_si(:), density_si(jtarg), density_si(jproj), regime, mlr, mslr, mtiny_si, Qloss)

            mass_res(1) = min(max(mlr, 0.0_DP), mtot)
            mass_res(2) = min(max(mslr, 0.0_DP), mtot)
            mass_res(3) = min(max(mtot - mlr - mslr, 0.0_DP), mtot)
            mass_res(:) = (mass_res(:) / MU2KG) * GU
            Qloss = Qloss * (GU / MU2KG) * (TU2S / DU2M)**2
         else !! When user has *not* enabled FRAGMENTATION, do every collision as a pure merger.
            regime = COLLRESOLVE_REGIME_MERGE
         end if

         write(*, *) "Collision detected at time t = ",t
         
         ! Set the appropriate flags for each of the discard types
         !! Use the positions and velocities of the parents and their children after the step is complete to generate the fragments
         select case (regime)
         case (COLLRESOLVE_REGIME_DISRUPTION)
            write(*, '("Disruption between particles ",20(I6,",",:))') parent_child_index_array(1)%id(:), parent_child_index_array(2)%id(:)
            status = symba_casedisruption(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss)
         case (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
            write(*, '("Supercatastrophic disruption between particles ",20(I6,",",:))') parent_child_index_array(1)%id(:), parent_child_index_array(2)%id(:)
            status = symba_casesupercatastrophic(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss)
         case (COLLRESOLVE_REGIME_HIT_AND_RUN)
            write(*, '("Hit and run between particles ",20(I6,",",:))') parent_child_index_array(1)%id(:), parent_child_index_array(2)%id(:)
            status = symba_casehitandrun(symba_plA, idx_parent, nmergeadd, mergeadd_list, id, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss)
         case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
            write(*, '("Merging particles ",20(I6,",",:))') parent_child_index_array(1)%id(:), parent_child_index_array(2)%id(:)
            status = symba_casemerge(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, param)
         case default 
            write(*,*) "Error in symba_collision, unrecognized collision regime"
            call util_exit(FAILURE)
            status = ACTIVE
         end select

         !write(*,*) 'Current status of all family members: '
         !write(*,*) ' index  id  status'
         !do i = 1, fam_size
         !   write(*,*) family(i),symba_plA%helio%swiftest%id(family(i)),symba_plA%helio%swiftest%status(family(i))
         !end do
         !write(*,*) 'Changing all status flags to: ',status

         ! If any body in the current collisional family is listed in subsequent collisions in this step, remove that 
         ! collision from consideration, as the body's outcome has already been resolved.
         where(statpl(family(:)) == COLLISION) statpl(family(:)) = status 
         do k = index_enc + 1, nplplenc
            if (plplenc_list%status(k) /= COLLISION) cycle ! Not the primary collision for this pair
            ! Index values of the original particle pair 
            idx(1) = plplenc_list%index1(k)
            idx(2) = plplenc_list%index2(k)
            if (any(family(:) == idx(1)) .or. any(family(:) == idx(2))) plplenc_list%status(k) = ACTIVE
         end do

        ! Msystem = sum(symba_plA%helio%swiftest%mass(1:npl), symba_plA%helio%swiftest%status(1:npl) == ACTIVE) 
        ! Mdiscard = sum(symba_plA%helio%swiftest%mass(1:npl), symba_plA%helio%swiftest%status(1:npl) /= ACTIVE) 
        ! Madd = sum(mergeadd_list%mass(1:nmergeadd))
        ! write(*,*) 'Mass balance in collision step'
        ! write(*,*) ' Msystem / Minitial:  ', Msystem / Minitial
        ! write(*,*) 'Mdiscard / Minitial: ', Mdiscard / Minitial
        ! write(*,*) '    Madd / Minitial: ', Madd / Minitial
        ! write(*,*) '  Mtotal / Minitial: ', (Msystem + Madd) / Minitial

         ! Reset the parent/child/family lists for the next collision
         do j = 1, 2
            deallocate(parent_child_index_array(j)%idx)
            deallocate(parent_child_index_array(j)%id)
         end do
         deallocate(family)

      end do
   end associate

   return

end subroutine symba_collision
