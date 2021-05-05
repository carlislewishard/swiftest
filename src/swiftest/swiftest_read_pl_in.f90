submodule (swiftest_data_structures) s_swiftest_read_pl_in
contains
   module procedure swiftest_read_pl_in !(self, param)
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      use swiftest
      use module_interfaces
      implicit none
      ! Arguments
      !class(swiftest_pl),          intent(inout) :: self   !! Swiftest data structure to store massive body initial condition
      !type(user_input_parameters), intent(inout) :: param  !! Input collection of user-defined parameters
      ! Internals
      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: i, ierr, npl
      logical                 :: is_ascii 

      ierr = 0
      is_ascii = (param%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = LUN, file = param%inplfile, status = 'old', form = 'formatted', iostat = ierr)
      else
         open(unit = LUN, file = param%inplfile, status = 'old', form = 'unformatted', iostat = ierr)
      end if
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(param%inplfile))
         return
      end if

      if (is_ascii) then
         read(LUN, *, iostat = ierr) npl
      else
         read(LUN, iostat = ierr) npl
      end if
      if (npl <= 0) return
      call self%alloc(npl)

      if (is_ascii) then
         read(LUN, *, iostat = ierr) self%id(1), self%mass(1)
         self%rhill(1) = 0.0_DP
         self%radius(1) = param%rmin
         read(LUN, *, iostat = ierr) self%xh(:,1)
         read(LUN, *, iostat = ierr) self%vh(:,1)
         if (param%lrotation) THEN
            read(LUN, *, iostat = ierr) self%Ip(:,1)
            read(LUN, *, iostat = ierr) self%rot(:,1)
         end if
         if (ierr /= 0) then
            write(*,*) 'Error reading central body values in ',trim(adjustl(param%inplfile))
            return
         end if
         do i = 1, NDIM
            if ((self%xh(i,1) /= 0.0_DP) .or. (self%vh(i,1) /= 0.0_DP)) then
               write(*, *) "Swiftest error:"
               write(*, *) " Input must be in heliocentric coordinates."
               write(*, *) " position/velocity components of body 1 are"
               write(*, *) self%xh(:,1)
               write(*, *) self%vh(:,1)
            end if
         end do
         self%status(1) = ACTIVE
         do i = 2, self%nbody
            if (param%lrhill_present) then
               read(LUN, *, iostat = ierr) self%id(i), self%mass(i), self%rhill(i)
            else
               read(LUN, *, iostat = ierr) self%id(i), self%mass(i)
               self%rhill(i) = 0.0_DP
            end if
            if (ierr /= 0 ) exit
            if (param%lclose) then
               read(LUN, *, iostat = ierr) self%radius(i)
               if (ierr /= 0 ) exit
            else
               self%radius(i) = 0.0_DP
            end if
            read(LUN, *, iostat = ierr) self%xh(:,i)
            read(LUN, *, iostat = ierr) self%vh(:,i)
            if (param%lrotation) THEN
               read(LUN, *, iostat = ierr) self%Ip(:,i)
               read(LUN, *, iostat = ierr) self%rot(:,i)
            end if
            if (ierr /= 0 ) exit
            self%status(i) = ACTIVE
         end do
      else
         read(LUN, iostat = ierr) self%id(:)
         read(LUN, iostat = ierr) self%mass(:)
         if (param%lrhill_present) then
            read(LUN, iostat = ierr) self%rhill(:)
         else
            self%rhill(:) = 0.0_DP
         end if
         self%status(:) = ACTIVE
         if (param%lclose) then
            read(LUN, iostat = ierr) self%radius(:)
         else
            self%radius(:) = 0.0_DP
         end if
         read(LUN, iostat = ierr) self%xh(:,:)
         read(LUN, iostat = ierr) self%vh(:,:)
         if (param%lrotation) THEN
            read(LUN, iostat = ierr) self%Ip(:,:)
            read(LUN, iostat = ierr) self%rot(:,:)
         end if
         self%status(:) = ACTIVE
      end if
      close(unit = LUN)
      if (ierr /= 0 ) then
         write(*,*) 'Error reading in massive body initial conditions from ',trim(adjustl(param%inplfile))
         call util_exit(FAILURE)
      end if

      do i = 1, npl
         self%info(i)%origin_xh(:) = self%xh(:,i)
         self%info(i)%origin_vh(:) = self%vh(:,i)
         self%info(i)%origin_time = 0.0_DP
         self%info(i)%origin_type = "Initial conditions"
      end do
      self%info(1)%origin_type = "Central body"

      ! Give massive bodies a positive id
      self%id(:) = abs(self%id(:))
      self%maxid = maxval(self%id(:))

      return
   end procedure swiftest_read_pl_in 

end submodule s_swiftest_read_pl_in

