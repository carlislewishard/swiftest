submodule (io) s_io_write_frame
contains
   module procedure io_write_frame
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Write a frame (header plus records for each plAnet and active test particle) to output binary file
   !! There is no direct file output from this subroutine
   !!
   !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
   !! Adapted from Hal Levison's Swift routine io_write_frame.F
   use module_interfaces
   implicit none
   

   logical                   :: lxdr
   logical, save             :: lfirst = .true.
   integer(I4B), parameter   :: lun = 20
   integer(I4B)              :: i, j, ierr
   integer(I4B), save        :: iu = lun, iout_form = xv
   real(DP)                  :: a, e, inc, capom, omega, capm, mu
   real(DP), dimension(NDIM) :: xtmp, vtmp
   real(DP), dimension(2:swiftest_plA%nbody)  :: a_pl, e_pl, inc_pl, capom_pl, omega_pl, capm_pl
   real(DP), dimension(1:swiftest_tpA%nbody)  :: a_tp, e_tp, inc_tp, capom_tp, omega_tp, capm_tp

   if (lfirst) then
      select case(out_stat)
      case('APPEND')
         open(unit = iu, file = outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      case('NEW')
         open(unit = iu, file = outfile, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      case ('REPLACE')
         open(unit = iu, file = outfile, status = 'REPLACE', form = 'UNFORMATTED', iostat = ierr)
      case default
         write(*,*) 'Invalid status code',trim(adjustl(out_stat))
         call util_exit(FAILURE)
      end select
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   binary output file already exists or cannot be accessed"
         call util_exit(FAILURE)
      end if

      select case (out_form)
      case ("EL")
         iout_form = EL
      case ("XV")
         iout_form = XV
      end select
      lfirst = .false.
   else
      open(unit = iu, file = outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   unable to open binary output file for APPEND"
         call util_exit(FAILURE)
      end if
   end if

   call io_write_hdr(iu, t, swiftest_plA%nbody, swiftest_tpA%nbody, iout_form, out_type)
   select case (iout_form)
   case (EL)
      do i = 2, swiftest_plA%nbody
         mu = swiftest_plA%mass(1) + swiftest_plA%mass(i)
         j = swiftest_plA%name(i)
         call orbel_xv2el(swiftest_plA%xh(:,i), swiftest_plA%vh(:,i), mu, a, e, inc, capom, omega, capm)
         a_pl(i) = a 
         e_pl(i) = e
         inc_pl(i) = inc
         capom_pl(i) = capom
         omega_pl(i) = omega
         capm_pl(i) = capm 
      end do
      mu = swiftest_plA%mass(1)
      do i = 1, swiftest_tpA%nbody
         j = swiftest_tpA%name(i)
         call orbel_xv2el(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, a, e, inc, capom, omega, capm)
         a_tp(i) = a 
         e_tp(i) = e
         inc_tp(i) = inc 
         capom_tp(i) = capom
         omega_tp(i) = omega
         capm_tp(i) = capm
      end do
      write(LUN) swiftest_plA%name(2:swiftest_plA%nbody)
      write(LUN) a_pl(2:swiftest_plA%nbody)
      write(LUN) e_pl(2:swiftest_plA%nbody)
      write(LUN) inc_pl(2:swiftest_plA%nbody)
      write(LUN) capom_pl(2:swiftest_plA%nbody)
      write(LUN) omega_pl(2:swiftest_plA%nbody)
      write(LUN) capm_pl(2:swiftest_plA%nbody)
      write(LUN) swiftest_plA%mass(2:swiftest_plA%nbody)
      write(LUN) swiftest_plA%radius(2:swiftest_plA%nbody)
      if (swiftest_tpA%nbody > 0) then
         write(LUN) swiftest_tpA%name(1:swiftest_tpA%nbody)
         write(LUN) a_tp(1:swiftest_tpA%nbody)
         write(LUN) e_tp(1:swiftest_tpA%nbody)
         write(LUN) inc_tp(1:swiftest_tpA%nbody)
         write(LUN) capom_tp(1:swiftest_tpA%nbody)
         write(LUN) omega_tp(1:swiftest_tpA%nbody)
         write(LUN) capm_tp(1:swiftest_tpA%nbody)
      end if

   case (XV)
         write(LUN) swiftest_pla%name(2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%xh(1,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%xh(2,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%xh(3,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%vh(1,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%vh(2,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%vh(3,2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%mass(2:swiftest_plA%nbody)
         write(LUN) swiftest_pla%radius(2:swiftest_plA%nbody) 

         if (swiftest_tpA%nbody > 0) then
            write(LUN) swiftest_tpa%name(:)  
            write(LUN) swiftest_tpa%xh(1,:)
            write(LUN) swiftest_tpa%xh(2,:)
            write(LUN) swiftest_tpa%xh(3,:)
            write(LUN) swiftest_tpa%vh(1,:)
            write(LUN) swiftest_tpa%vh(2,:)
            write(LUN) swiftest_tpa%vh(3,:)
         end if 

   end select


   close(unit = iu, iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   unable to close binary output file"
      call util_exit(FAILURE)
   end if

   return

   end procedure io_write_frame
end submodule s_io_write_frame
