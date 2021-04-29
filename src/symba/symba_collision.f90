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
   integer(I4B)                            :: i, j, index_enc, jtarg, jproj
   real(DP), dimension(NRES)               :: mass_res
   real(DP), dimension(NDIM)               :: x1_si, v1_si, x2_si, v2_si
   integer(I4B)                            :: regime, idx_child, status
   integer(I4B), dimension(2)              :: idx, idx_parent, nchild, name 
   real(DP), dimension(2)                  :: radius, mass, density, volume
   real(DP), dimension(2)                  :: radius_si, mass_si, density_si
   real(DP), dimension(NDIM, 2)            :: x, v, L_spin, Ip
   real(DP)                                :: volchild, dentot, Mcb_si
   real(DP)                                :: mchild, mtot
   real(DP), dimension(NDIM)               :: xc, vc, xcom, vcom, xchild, vchild, xcrossv
   real(DP)                                :: mtiny_si
   integer(I4B), dimension(:), allocatable :: array_index1_child, array_index2_child, name1, name2
   real(DP)                                :: mlr, mslr, msys, msys_new, Qloss
   logical                                 :: lpure


   ! First determine the collisional regime for each colliding pair
   associate(npl => symba_plA%helio%swiftest%nbody, xbpl => symba_plA%helio%swiftest%xb)
      ldiscard = any(plplenc_list%status(1:nplplenc) == COLLISION)
      if (.not.ldiscard) return

      ! Recompute central body barycentric velocity
      call coord_h2b(npl, symba_plA%helio%swiftest, msys)

      ! Set the appropriate flags for each of the discard types
      do index_enc = 1, nplplenc
         if (plplenc_list%status(index_enc) /= COLLISION) cycle ! Not the primary collision for this pair

         ! Index values of the original particle pair 
         idx(1) = plplenc_list%index1(index_enc)
         idx(2) = plplenc_list%index2(index_enc)

         if (any(symba_plA%helio%swiftest%status(idx(:)) /= ACTIVE)) cycle ! One of these two bodies is already gone
         symba_plA%helio%swiftest%status(idx(:)) = COLLISION

         ! Index values for the parents of this particle pair
         idx_parent(:) = symba_plA%kin(idx(:))%parent

         nchild(:) = symba_plA%kin(idx_parent(:))%nchild 
         ! If all of these bodies share a parent, but this is still a unique collision, move the last child
         ! out of the parent's position and make it the secondary body
         if (idx_parent(1) == idx_parent(2)) then
            idx_parent(2) = symba_plA%kin(idx_parent(1))%child(nchild(1))
            nchild(1) = nchild(1) - 1
            symba_plA%kin(idx_parent(1))%nchild = nchild(1)
            nchild(2) = 0
         end if

         mass(:) = symba_plA%helio%swiftest%mass(idx_parent(:))
         name(:) = symba_plA%helio%swiftest%name(idx_parent(:))
         radius(:) = symba_plA%helio%swiftest%radius(idx_parent(:))
         volume(:) =  (4.0_DP / 3.0_DP) * PI * radius(:)**3
      
         if (nchild(1) > 0) then
            allocate(array_index1_child, source = symba_plA%kin(idx_parent(1))%child(1:nchild(1)))
            allocate(name1(nchild(1)+1))
         else 
            allocate(array_index1_child(1))
            allocate(name1(1))
            array_index1_child(1) = idx_parent(1) 
         end if
      
         if (nchild(2) > 0) then
            allocate(array_index2_child, source = symba_plA%kin(idx_parent(2))%child(1:nchild(2)))
            allocate(name2(nchild(2)+1))
         else 
            allocate(array_index2_child(1))
            allocate(name2(1))
            array_index2_child(1) = idx_parent(2)
         end if
         name1(1) = name(1)
         name2(1) = name(2)

         ! Find the barycenter of each body along with its children, if it has any
         do j = 1, 2
            x(:, j)  = symba_plA%helio%swiftest%xb(:, idx_parent(j))
            v(:, j)  = symba_plA%helio%swiftest%vb(:, idx_parent(j))
            Ip(:, j) = mass(j) * symba_plA%helio%swiftest%Ip(:, idx_parent(j))
            ! Assume principal axis rotation about axis corresponding to highest moment of inertia (3rd Ip)
            L_spin(:, j)  = Ip(3, j) * radius(j)**2 * symba_plA%helio%swiftest%rot(:, idx_parent(j))
            if (nchild(j) > 0) then
               do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
                  if (j == 1) then
                     idx_child = array_index1_child(i)
                     name1(1+i) = symba_plA%helio%swiftest%name(idx_child)
                  else
                     idx_child = array_index2_child(i)
                     name2(1+i) = symba_plA%helio%swiftest%name(idx_child)
                  end if
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
         status = COLLISION
         ! Set the appropriate flags for each of the discard types
         !! Use the positions and velocities of the parents and their children after the step is complete to generate the fragments
         select case (regime)
         case (COLLRESOLVE_REGIME_DISRUPTION)
            write(*, '("Disruption between particles ",20(I6,",",:))') name1(:), name2(:) 
            call symba_casedisruption(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss)
         case (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
            write(*, '("Supercatastrophic disruption between particles ",20(I6,",",:))') name1(:), name2(:) 
            call symba_casesupercatastrophic(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss)
         case (COLLRESOLVE_REGIME_HIT_AND_RUN)
            write(*, '("Hit and run between particles ",20(I6,",",:))') name1(:), name2(:) 
            call symba_casehitandrun(symba_plA, idx_parent, nmergeadd, mergeadd_list, name, x, v, mass, radius, L_spin, Ip, mass_res, param, Qloss, lpure)
            if (lpure) status = ACTIVE
         case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
            write(*, '("Merging particles ",20(I6,",",:))') name1(:), name2(:) 
            call symba_casemerge(symba_plA, idx_parent, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, param)
         case default 
            write(*,*) "Error in symba_collision, unrecognized collision regime"
            call util_exit(FAILURE)
         end select

         symba_plA%helio%swiftest%status(idx_parent(:)) = status
         do j = 1, 2
            if (nchild(j) > 0) then
               do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
                  if (j == 1) then
                     idx_child = array_index1_child(i)
                  else
                     idx_child = array_index2_child(i)
                  end if
                  symba_plA%helio%swiftest%status(idx_child) = status 
               end do
            end if
         end do
         deallocate(array_index1_child, array_index2_child, name1, name2)
      end do
   end associate

   return

end subroutine symba_collision
