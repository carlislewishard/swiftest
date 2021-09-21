submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains

   subroutine symba_encounter_check_all(nplplm, k_plpl, x, v, rhill,  dt, irec, lencounter, loc_lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance
      implicit none
      integer(I8B), intent(in) :: nplplm
      integer(I4B), dimension(:,:), intent(in) :: k_plpl
      real(DP), dimension(:,:), intent(in) :: x, v
      real(DP), dimension(:), intent(in) :: rhill
      real(DP), intent(in) :: dt
      integer(I4B), intent(in) :: irec
      logical, dimension(:), intent(out) :: lencounter, loc_lvdotr
      ! Internals
      integer(I8B) :: k
      integer(I4B) :: i, j
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2

      !$omp parallel do default(private)&
      !$omp shared(nplplm, k_plpl, x, v, rhill, dt, irec, lencounter, loc_lvdotr)
      do k = 1_I8B, nplplm
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         xr = x(1, j) - x(1, i)
         yr = x(2, j) - x(2, i)
         zr = x(3, j) - x(3, i)
         vxr = v(1, j) - v(1, i)
         vyr = v(2, j) - v(2, i)
         vzr = v(3, j) - v(3, i)
         rhill1 = rhill(i)
         rhill2 = rhill(j)
         lencounter(k) = .false.
         loc_lvdotr(k) = .false.
         call symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter(k), loc_lvdotr(k))
      end do
      !$omp end parallel do

      return
   end subroutine symba_encounter_check_all


   module function symba_encounter_check_pl(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout)  :: self           !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout)  :: system         !! SyMBA nbody system object
      real(DP),                  intent(in)     :: dt             !! step size
      integer(I4B),              intent(in)     :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I8B) :: k, nplplm
      integer(I4B) :: i, j, nenc
      logical, dimension(:), allocatable :: lencounter, loc_lvdotr
  
      if (self%nbody == 0) return

      associate(pl => self)
         nplplm = pl%nplplm
         allocate(lencounter(nplplm))
         allocate(loc_lvdotr(nplplm))
  
         call symba_encounter_check_all(nplplm, pl%k_plpl, pl%xh, pl%vh, pl%rhill, dt, irec, lencounter, loc_lvdotr)

         !$omp parallel workshare
         nenc = count(lencounter(:))
         !$omp end parallel workshare

         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(plplenc_list => system%plplenc_list)
               call plplenc_list%resize(nenc)
               plplenc_list%lvdotr(1:nenc) = pack(loc_lvdotr(1:nplplm), lencounter(1:nplplm))
               plplenc_list%kidx(1:nenc) = pack([(k, k = 1_I8B, nplplm)], lencounter(1:nplplm))
               deallocate(lencounter, loc_lvdotr)
               plplenc_list%index1(1:nenc) = pl%k_plpl(1,plplenc_list%kidx(1:nenc))
               plplenc_list%index2(1:nenc) = pl%k_plpl(2,plplenc_list%kidx(1:nenc))
               plplenc_list%id1(1:nenc) = pl%id(plplenc_list%index1(1:nenc))
               plplenc_list%id2(1:nenc) = pl%id(plplenc_list%index2(1:nenc))
               do k = 1, nenc
                  i = plplenc_list%index1(k)
                  j = plplenc_list%index2(k)
                  plplenc_list%status(k) = ACTIVE
                  plplenc_list%level(k) = irec
                  pl%lencounter(i) = .true.
                  pl%lencounter(j) = .true.
                  pl%levelg(i) = irec
                  pl%levelm(i) = irec
                  pl%levelg(j) = irec
                  pl%levelm(j) = irec
                  pl%nplenc(i) = pl%nplenc(i) + 1
                  pl%nplenc(j) = pl%nplenc(j) + 1
               end do
            end associate
         end if
      end associate
      return
   end function symba_encounter_check_pl


   module function symba_encounter_check(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies in the pltpenc list.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types.
      !!
      !! Adapted from portions of David E. Kaufmann's Swifter routine: symba_step_recur.f90
      implicit none
      ! Arguments
      class(symba_encounter),      intent(inout) :: self       !! SyMBA pl-pl encounter list object
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec       !! Current recursion level 
      logical                                  :: lany_encounter !! Returns true if there is at least one close encounter  
      ! Internals
      integer(I4B)              :: i, j,k
      real(DP), dimension(NDIM) :: xr, vr
      logical                   :: lencounter, isplpl
      real(DP)                  :: rlim2, rji2
      logical, dimension(:), allocatable :: lencmask

      lany_encounter = .false.
      if (self%nenc == 0) return

      select type(self)
      class is (symba_plplenc)
         isplpl = .true.
      class is (symba_pltpenc)
         isplpl = .false.
      end select

      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            allocate(lencmask(self%nenc))
            lencmask(:) = (self%status(1:self%nenc) == ACTIVE) .and. (self%level(1:self%nenc) == irec - 1)
            if (.not.any(lencmask(:))) return
            do concurrent(k = 1:self%nenc, lencmask(k))
               i = self%index1(k)
               j = self%index2(k)
               if (isplpl) then
                  xr(:) = pl%xh(:,j) - pl%xh(:,i)
                  vr(:) = pl%vb(:,j) - pl%vb(:,i)
                  call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), pl%rhill(j), dt, irec, lencounter, self%lvdotr(k))
               else
                  xr(:) = tp%xh(:,j) - pl%xh(:,i)
                  vr(:) = tp%vb(:,j) - pl%vb(:,i)
                  call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), 0.0_DP, dt, irec, lencounter, self%lvdotr(k))
               end if
               if (lencounter) then
                  if (isplpl) then
                     rlim2 = (pl%radius(i) + pl%radius(j))**2
                  else
                     rlim2 = (pl%radius(i))**2
                  end if
                  rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                  if (rji2 > rlim2) then
                     lany_encounter = .true.
                     pl%levelg(i) = irec
                     pl%levelm(i) = MAX(irec, pl%levelm(i))
                     if (isplpl) then
                        pl%levelg(j) = irec
                        pl%levelm(j) = MAX(irec, pl%levelm(j))
                     else
                        tp%levelg(j) = irec
                        tp%levelm(j) = MAX(irec, tp%levelm(j))
                     end if
                     self%level(k) = irec
                  end if
               end if   
            end do
         end select
      end select

      return
   end function symba_encounter_check


   module function symba_encounter_check_tp(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      real(DP)                                  :: r2crit, vdotr, r2, v2, tmin, r2min, term2
      integer(I4B)                              :: i, j, k,nenc, plind, tpind
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:,:), allocatable :: lencounter, loc_lvdotr
  
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => system%pl, npl => system%pl%nbody)
         allocate(lencounter(ntp, npl), loc_lvdotr(ntp, npl))
         lencounter(:,:) = .false.
   
         do j = 1, npl
            do i = 1, ntp
               xr(:) = tp%xh(:, i) - pl%xh(:, j)
               vr(:) = tp%vh(:, i) - pl%vh(:, j)
               call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(j), 0.0_DP, dt, irec, lencounter(i,j), loc_lvdotr(i,j))
            end do
         end do

         nenc = count(lencounter(:,:))
         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(pltpenc_list => system%pltpenc_list)
               call pltpenc_list%resize(nenc)
               pltpenc_list%status(1:nenc) = ACTIVE
               pltpenc_list%level(1:nenc) = irec
               pltpenc_list%lvdotr(1:nenc) = pack(loc_lvdotr(1:ntp, 1:npl), lencounter(1:ntp, 1:npl))
               pltpenc_list%index1(1:nenc) = pack(spread([(i, i = 1, npl)], dim=1, ncopies=ntp), lencounter(1:ntp, 1:npl)) 
               pltpenc_list%index2(1:nenc) = pack(spread([(i, i = 1, ntp)], dim=2, ncopies=npl), lencounter(1:ntp, 1:npl))
               pltpenc_list%id1(1:nenc) = pl%id(pltpenc_list%index1(1:nenc))
               pltpenc_list%id2(1:nenc) = tp%id(pltpenc_list%index2(1:nenc))
               select type(pl)
               class is (symba_pl)
                  pl%lencounter(1:npl) = .false.
                  do k = 1, nenc
                     plind = pltpenc_list%index1(k)
                     tpind = pltpenc_list%index2(k)
                     pl%lencounter(plind) = .true.
                     pl%levelg(plind) = irec
                     pl%levelm(plind) = irec
                     tp%levelg(tpind) = irec
                     tp%levelm(tpind) = irec
                     pl%ntpenc(plind) = pl%ntpenc(plind) + 1
                     tp%nplenc(tpind) = tp%nplenc(tpind) + 1
                  end do
               end select
            end associate
         end if
      end associate

      return
   end function symba_encounter_check_tp


   module pure elemental subroutine symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for an encounter.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_chk.f90
      !! Adapted from Hal Levison's Swift routine symba5_chk.f
      implicit none
      ! Arguments
      real(DP),     intent(in)  :: xr, yr, zr     !! Relative distance vector components
      real(DP),     intent(in)  :: vxr, vyr, vzr  !! Relative velocity vector components
      real(DP),     intent(in)  :: rhill1, rhill2 !! Hill spheres of the two bodies
      real(DP),     intent(in)  :: dt             !! Step size
      integer(I4B), intent(in)  :: irec           !! Current SyMBA recursion level
      real(DP),     intent(in)  :: r2crit         !! Square of the critical encounter distance
      logical,      intent(out) :: lencounter     !! Flag indicating that an encounter has occurred
      logical,      intent(out) :: lvdotr         !! Logical flag indicating the direction of the v .dot. r vector
      ! Internals
      real(DP)     :: r2crit

      r2crit = (rhill1 + rhill2)*RHSCALE*(RSHELL**(irec))
      r2crit = r2crit**2
      call rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, lencounter, lvdotr)

      return
   end subroutine symba_encounter_check_one

end submodule s_symba_encounter_check