subroutine io_dump_pl_symba(npl, symba_plA, param)
   !! Author: David A. Minton
   !!
   !! Dumps a SyMBA data structure to file.
   !! 
! modules
   use swiftest
   use module_symba
   use module_interfaces, except_this_one => io_dump_pl_symba
   implicit none

! arguments
   integer(I4B), intent(in)         :: npl
   type(symba_pl), intent(inout) :: symba_plA
   type(io_input_parameters),intent(inout) :: param

! internals
   integer(I4B)                     :: ierr
   integer(I4B), save               :: idx = 1
   integer(I4B),parameter           :: lun = 7

   open(unit = lun, file = dump_pl_file(idx), form = "UNFORMATTED", status = 'REPLACE', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Unable to open binary dump file ", trim(dump_pl_file(idx))
      call util_exit(FAILURE)
   end if
   write(lun) npl
   write(lun) symba_plA%helio%swiftest%id(1:npl)
   write(lun) symba_plA%helio%swiftest%mass(1:npl)
   if (param%lrhill_present) write(lun) symba_plA%helio%swiftest%rhill(1:npl) 
   if (param%lclose) write(lun) symba_plA%helio%swiftest%radius(1:npl) 
   write(lun) symba_plA%helio%swiftest%xh(:,1:npl)
   write(lun) symba_plA%helio%swiftest%vh(:,1:npl)
   if (param%lrotation) then
      write(lun) symba_plA%helio%swiftest%ip(:,1:npl)
      write(lun) symba_plA%helio%swiftest%rot(:,1:npl)
   end if
   write(lun) symba_plA%helio%swiftest%vb(:,1:npl)
   write(lun) symba_plA%helio%ah(:,1:npl)
   write(lun) symba_plA%helio%ahi(:,1:npl)
   close(lun)
   if (idx == 1) then
      idx = 2
   else
      idx = 1
   end if

   return

end subroutine io_dump_pl_symba