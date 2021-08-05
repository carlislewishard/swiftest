submodule (symba_classes) s_symba_fragmentation
   use swiftest

   integer(I4B), parameter :: NFRAG_DISRUPT = 12
   integer(I4B), parameter :: NFRAG_SUPERCAT = 20
contains

   module function symba_fragmentation_casedisruption(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, istart, nfrag, ibiggest, nstart, nend
      real(DP)                                :: mtot, avg_dens
      real(DP), dimension(NDIM)               :: xcom, vcom, Ip_new
      real(DP), dimension(2)                  :: vol
      real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
      real(DP), dimension(:), allocatable     :: m_frag, rad_frag
      logical                                 :: lfailure
      class(symba_pl), allocatable            :: plnew
   
      select type(pl => system%pl)
      class is (symba_pl)
         associate(mergeadd_list => system%mergeadd_list, mergesub_list => system%mergesub_list, cb => system%cb)
            ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
            nfrag = NFRAG_DISRUPT 
            allocate(m_frag(nfrag))
            allocate(rad_frag(nfrag))
            allocate(xb_frag(NDIM, nfrag))
            allocate(vb_frag(NDIM, nfrag))
            allocate(rot_frag(NDIM, nfrag))
            allocate(Ip_frag(NDIM, nfrag))
         
            mtot = sum(mass(:))
            xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
            vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot
         
            ! Get mass weighted mean of Ip and average density
            Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mtot
            vol(:) = 4._DP / 3._DP * PI * radius(:)**3
            avg_dens = mtot / sum(vol(:))
         
            ! Distribute the mass among fragments, with a branch to check for the size of the second largest fragment
            m_frag(1) = mass_res(1)
            if (mass_res(2) > mass_res(1) / 3._DP) then
               m_frag(2) = mass_res(2)
               istart = 3
            else
               istart = 2
            end if
            ! Distribute remaining mass among the remaining bodies
            do i = istart, nfrag
               m_frag(i) = (mtot - sum(m_frag(1:istart - 1))) / (nfrag - istart + 1) 
            end do
         
            ! Distribute any residual mass if there is any and set the radius
            m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
            rad_frag(:) = (3 * m_frag(:) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)
         
            do i = 1, nfrag
               Ip_frag(:, i) = Ip_new(:)
            end do
         
            !call fragmentation_initialize(pl, param, family, x, v, L_spin, Ip, mass, radius, &
            !                              nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lfailure)
         
            if (lfailure) then
               write(*,*) 'No fragment solution found, so treat as a pure hit-and-run'
               status = ACTIVE 
               nfrag = 0
            else
               ! Populate the list of new bodies
               write(*,'("Generating ",I2.0," fragments")') nfrag
               allocate(plnew, mold=pl)
               call plnew%setup(nfrag, param)
    
               plnew%id(:) = [(i, i = system%maxid + 1, system%maxid + nfrag)]
               system%maxid = system%maxid + nfrag
               plnew%status(:) = ACTIVE
               plnew%lcollision(:) = .false.
               plnew%ldiscard(:) = .false.
               plnew%xb(:,:) = xb_frag(:, :) 
               plnew%vb(:,:) = vb_frag(:, :)
               do i = 1, nfrag
                  plnew%xh(:,i) = xb_frag(:, i) - cb%xb(:)
                  plnew%vh(:,i) = vb_frag(:, i) - cb%vb(:)
               end do
               plnew%mass(:) = m_frag(:)
               plnew%Gmass(:) = param%GU * m_frag(:)
               plnew%density(:) = avg_dens
               plnew%radius(:) = rad_frag(:)
               plnew%info(:)%origin_type = "Disruption"
               plnew%info(:)%origin_time = param%t
               do i = 1, nfrag
                  plnew%info(i)%origin_xh(:) = plnew%xh(:,i)
                  plnew%info(i)%origin_vh(:) = plnew%vh(:,i)
               end do
               if (param%lrotation) then
                  plnew%Ip(:,:) = Ip_frag(:,:)
                  plnew%rot(:,:) = rot_frag(:,:)
               end if
               if (param%ltides) then
                  ibiggest = maxloc(pl%Gmass(family(:)), dim=1)
                  plnew%Q = pl%Q(ibiggest)
                  plnew%k2 = pl%k2(ibiggest)
                  plnew%tlag = pl%tlag(ibiggest)
               end if
   
               ! Append the new merged body to the list and record how many we made
               nstart = mergeadd_list%nbody + 1
               nend = mergeadd_list%nbody + plnew%nbody
               call mergeadd_list%append(plnew)
               mergeadd_list%ncomp(nstart:nend) = plnew%nbody
    
               call plnew%setup(0, param)
               deallocate(plnew)
            end if

         end associate
      end select

      return
   end function symba_fragmentation_casedisruption


   module function symba_fragmentation_casehitandrun(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals

      status = HIT_AND_RUN

      return
   end function symba_fragmentation_casehitandrun


   module function symba_fragmentation_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge planets.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(in)    :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(in)    :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                              :: i, j, ibiggest, nfamily, nstart, nend
      real(DP)                                  :: mass_new, radius_new, volume_new, pe
      real(DP), dimension(NDIM)                 :: xcom, vcom, xc, vc, xcrossv
      real(DP), dimension(2)                    :: vol
      real(DP), dimension(NDIM)                 :: L_orb_old, L_spin_old
      real(DP), dimension(NDIM)                 :: L_spin_new, rot_new, Ip_new
      logical,  dimension(system%pl%nbody)      :: lmask
      class(symba_pl), allocatable              :: plnew
   
      select type(pl => system%pl)
      class is (symba_pl)
         associate(mergeadd_list => system%mergeadd_list, mergesub_list => system%mergesub_list, cb => system%cb)
            status = MERGED
            write(*, '("Merging bodies ",99(I8,",",:))') pl%id(family(:))
            mass_new = sum(mass(:))
      
            ! Merged body is created at the barycenter of the original bodies
            xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mass_new
            vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mass_new
      
            ! Get mass weighted mean of Ip and 
            vol(:) = 4._DP / 3._DP * PI * radius(:)**3
            volume_new = sum(vol(:))
            radius_new = (3 * volume_new / (4 * PI))**(1._DP / 3._DP)

            L_orb_old(:) = 0.0_DP

            ! Compute orbital angular momentum of pre-impact system
            do i = 1, 2
               xc(:) = x(:, i) - xcom(:)
               vc(:) = v(:, i) - vcom(:)
               xcrossv(:) = xc(:) .cross. vc(:)
               L_orb_old(:) = L_orb_old(:) + mass(i) * xcrossv(:)
            end do
         
            if (param%lrotation) then
               Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mass_new
               L_spin_old(:) = L_spin(:,1) + L_spin(:,2)

               ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
               L_spin_new(:) = L_orb_old(:) + L_spin_old(:) 
      
               ! Assume prinicpal axis rotation on 3rd Ip axis
               rot_new(:) = L_spin_new(:) / (Ip_new(3) * mass_new * radius_new**2)
            else ! If spin is not enabled, we will consider the lost pre-collision angular momentum as "escaped" and add it to our bookkeeping variable
               system%Lescape(:) = system%Lescape(:) + L_orb_old(:) 
            end if
      
            ! Keep track of the component of potential energy due to the pre-impact family for book-keeping
            nfamily = size(family(:))
            pe = 0.0_DP
            do j = 1, nfamily
               do i = j + 1, nfamily
                  pe = pe - pl%mass(i) * pl%mass(j) / norm2(pl%xb(:, i) - pl%xb(:, j))
               end do
            end do
            system%Ecollisions  = system%Ecollisions + pe 
            system%Euntracked = system%Euntracked - pe 
     
            ! Add the family bodies to the subtraction list
            lmask(:) = .false.
            lmask(family(:)) = .true.
            pl%status(family(:)) = MERGED
            nstart = mergesub_list%nbody + 1
            nend = mergesub_list%nbody + nfamily
            call mergesub_list%append(pl, lmask)
            ! Record how many bodies were subtracted in this event
            mergesub_list%ncomp(nstart:nend) = nfamily 

            ! Create the new merged body 
            allocate(plnew, mold=pl)
            call plnew%setup(1, param)

            ! The merged body's name will be that of the largest of the two parents 
            ibiggest = maxloc(pl%Gmass(family(:)), dim=1)
            plnew%id(1) = pl%id(family(ibiggest))
            plnew%status(1) = ACTIVE
            plnew%lcollision = .false.
            plnew%ldiscard = .false.
            plnew%xb(:,1) = xcom(:)
            plnew%vb(:,1) = vcom(:)
            plnew%xh(:,1) = xcom(:) - cb%xb(:)
            plnew%vh(:,1) = vcom(:) - cb%vb(:)
            plnew%mass(1) = mass_new
            plnew%Gmass(1) = param%GU * mass_new
            plnew%density(1) = mass_new / volume_new
            plnew%radius(1) = radius_new
            plnew%info(1) = pl%info(family(ibiggest)) 
            if (param%lrotation) then
               plnew%Ip(:,1) = Ip_new(:)
               plnew%rot(:,1) = rot_new(:)
            end if
            if (param%ltides) then
               plnew%Q = pl%Q(ibiggest)
               plnew%k2 = pl%k2(ibiggest)
               plnew%tlag = pl%tlag(ibiggest)
            end if

            ! Append the new merged body to the list and record how many we made
            nstart = mergeadd_list%nbody + 1
            nend = mergeadd_list%nbody + plnew%nbody
            call mergeadd_list%append(plnew)
            mergeadd_list%ncomp(nstart:nend) = plnew%nbody

            call plnew%setup(0, param)
            deallocate(plnew)

         end associate
      end select
   
      return 

   end function symba_fragmentation_casemerge


   module function symba_fragmentation_casesupercatastrophic(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a supercatastrophic collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals

      status = SUPERCATASTROPHIC

      return
   end function symba_fragmentation_casesupercatastrophic

end submodule s_symba_fragmentation