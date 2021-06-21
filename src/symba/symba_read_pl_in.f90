submodule (module_symba) s_symba_read_pl_in
contains
   module subroutine symba_read_pl_in(symba_plA, param) 
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine symba_init_pl.f
      use swiftest
      use module_interfaces
      implicit none
      ! Arguments
      type(symba_pl),              intent(inout) :: symba_plA  !! Swiftest data structure to store massive body initial conditions
      type(user_input_parameters), intent(inout) :: param    !! Input collection of user-defined parameters
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
      call symba_plA%helio%swiftest%alloc(npl)

      if (is_ascii) then
         read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%id(1), symba_plA%helio%swiftest%mass(1)
         symba_plA%helio%swiftest%rhill(1) = 0.0_DP
         symba_plA%helio%swiftest%radius(1) = param%rmin
         read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%xh(:,1)
         read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%vh(:,1)
         if (param%lrotation) THEN
            read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%Ip(:,1)
            read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%rot(:,1)
         end if
         if (ierr /= 0) then
            write(*,*) 'Error reading central body values in ',trim(adjustl(param%inplfile))
            return
         end if
         do i = 1, NDIM
            if ((symba_plA%helio%swiftest%xh(i,1) /= 0.0_DP) .or. (symba_plA%helio%swiftest%vh(i,1) /= 0.0_DP)) then
               write(*, *) "Swiftest error:"
               write(*, *) " Input must be in heliocentric coordinates."
               write(*, *) " position/velocity components of body 1 are"
               write(*, *) symba_plA%helio%swiftest%xh(:,1)
               write(*, *) symba_plA%helio%swiftest%vh(:,1)
            end if
         end do
         symba_plA%helio%swiftest%status(1) = ACTIVE
         do i = 2, symba_plA%helio%swiftest%nbody
            if (param%lrhill_present) then
               read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%id(i), symba_plA%helio%swiftest%mass(i), symba_plA%helio%swiftest%rhill(i)
            else
               read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%id(i), symba_plA%helio%swiftest%mass(i)
               symba_plA%helio%swiftest%rhill(i) = 0.0_DP
            end if
            if (ierr /= 0 ) exit
            if (param%lclose) then
               read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%radius(i)
               if (ierr /= 0 ) exit
            else
               symba_plA%helio%swiftest%radius(i) = 0.0_DP
            end if
            read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%xh(:,i)
            read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%vh(:,i)
            if (param%lrotation) THEN
               read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%Ip(:,i)
               read(LUN, *, iostat = ierr) symba_plA%helio%swiftest%rot(:,i)
            end if
            if (ierr /= 0 ) exit
            symba_plA%helio%swiftest%status(i) = ACTIVE
         end do
      else
         read(LUN, iostat = ierr) symba_plA%helio%swiftest%id(:)
         read(LUN, iostat = ierr) symba_plA%helio%swiftest%mass(:)
         if (param%lrhill_present) then
            read(LUN, iostat = ierr) symba_plA%helio%swiftest%rhill(:)
         else
            symba_plA%helio%swiftest%rhill(:) = 0.0_DP
         end if
         if (param%lclose) then
            read(LUN, iostat = ierr) symba_plA%helio%swiftest%radius(:)
         else
            symba_plA%helio%swiftest%radius(:) = 0.0_DP
         end if
         read(LUN, iostat = ierr) symba_plA%helio%swiftest%xh(:,:)
         read(LUN, iostat = ierr) symba_plA%helio%swiftest%vh(:,:)
         if (param%lrotation) THEN
            read(LUN, iostat = ierr) symba_plA%helio%swiftest%Ip(:,:)
            read(LUN, iostat = ierr) symba_plA%helio%swiftest%rot(:,:)
         end if
         if (param%out_stat == 'APPEND') then
            read(lun, iostat = ierr) symba_plA%helio%swiftest%vb(:,:)
            read(lun, iostat = ierr) symba_plA%helio%ah(:,:)
            read(lun, iostat = ierr) symba_plA%helio%ahi(:,:)
         end if
         if (ierr /= 0) then
            write(*,*) 'An error occurred reading in ',trim(adjustl(param%inplfile))
            call util_exit(FAILURE)
         end if
         symba_plA%helio%swiftest%status(:) = ACTIVE
      end if
      close(unit = LUN)

      symba_plA%helio%swiftest%info(1)%origin_type = "Central body"
      do i = 2, npl
         symba_plA%helio%swiftest%info(i)%origin_xh(:) = symba_plA%helio%swiftest%xh(:,i)
         symba_plA%helio%swiftest%info(i)%origin_vh(:) = symba_plA%helio%swiftest%vh(:,i)
         symba_plA%helio%swiftest%info(i)%origin_time = 0.0_DP
         symba_plA%helio%swiftest%info(i)%origin_type = "Initial conditions"
      end do

      ! Give massive bodies a positive id
      symba_plA%helio%swiftest%id(:) = abs(symba_plA%helio%swiftest%id(:))
      symba_plA%helio%swiftest%maxid = maxval(symba_plA%helio%swiftest%id(:))

      return
   end subroutine symba_read_pl_in

end submodule s_symba_read_pl_in

