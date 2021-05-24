submodule (io) s_io_read_particle

contains
   module subroutine io_read_particle_pl(swiftest_plA, param)
      !! author: David A. Minton
      !!
      !! Writes particle information to a file.
      !!
      use swiftest, except_this_one => io_read_particle_pl
      implicit none
      ! Arguments
      class(swiftest_pl),          intent(inout) :: swiftest_plA   !! Swiftest massive body structure
      type(user_input_parameters), intent(in) :: param   !! Input colleciton of user-defined parameters
      ! Internals
      integer(I4B), parameter   :: lun = 22
      integer(I4B)              :: i, ierr, id, idx
      integer(I4B), save        :: iu = lun

      associate(out_form => param%out_form, out_type => param%out_type, particle_file => param%particle_file)
         open(unit = iu, file = particle_file, status = 'OLD', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   unable to open binary particle file for reading"
            call util_exit(FAILURE)
         end if

         do 
            read(LUN, iostat=ierr) id
            if (ierr /=0) exit
            idx = findloc(swiftest_plA%id(:), id, dim=1)
            read(LUN) swiftest_plA%info(idx)
         end do
         close(unit = iu, iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   unable to close particle output file"
            call util_exit(FAILURE)
         end if
      end associate
      return

      end subroutine io_read_particle_pl
end submodule s_io_read_particle
