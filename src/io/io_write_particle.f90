submodule (io) s_io_write_particle

contains
   module subroutine io_write_particle_pl(swiftest_plA, idx, param)
      !! author: David A. Minton
      !!
      !! Writes particle information to a file.
      !!
      use swiftest, except_this_one => io_write_particle_pl
      implicit none
      ! Arguments
      class(swiftest_pl),          intent(in) :: swiftest_plA !! Swiftest massive body structure
      integer(I4B), dimension(:),  intent(in) :: idx       !! Array of particle indices to append to the particle file
      type(user_input_parameters), intent(in) :: param     !! Input colleciton of user-defined parameters
      ! Internals
      logical, save             :: lfirst = .true.
      integer(I4B), parameter   :: lun = 22
      integer(I4B)              :: i, ierr
      integer(I4B), save        :: iu = lun

      associate(out_stat => param%out_stat, out_form => param%out_form, out_type => param%out_type, particle_file => param%particle_file)
         if (lfirst) then
            select case(out_stat)
            case('APPEND')
               open(unit = iu, file = particle_file, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
            case('NEW')
               open(unit = iu, file = particle_file, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
            case ('REPLACE')
               open(unit = iu, file = particle_file, status = 'REPLACE', form = 'UNFORMATTED', iostat = ierr)
            case default
               write(*,*) 'Invalid status code',trim(adjustl(out_stat))
               call util_exit(FAILURE)
            end select
            if (ierr /= 0) then
               write(*, *) "Swiftest error:"
               write(*, *) "   particle output file already exists or cannot be accessed"
               call util_exit(FAILURE)
            end if

            lfirst = .false.
         else
            open(unit = iu, file = particle_file, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
            if (ierr /= 0) then
               write(*, *) "Swiftest error:"
               write(*, *) "   unable to open binary output file for APPEND"
               call util_exit(FAILURE)
            end if
         end if

         do i = 1, size(idx)
            write(LUN) swiftest_plA%id(idx(i))
            write(LUN) swiftest_plA%info(idx(i))
         end do
         close(unit = iu, iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   unable to close particle output file"
            call util_exit(FAILURE)
         end if
      end associate
      return

      end subroutine io_write_particle_pl
end submodule s_io_write_particle
