submodule(user) s_user_dump_param
contains
   module procedure user_dump_param
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      use swiftest
      use module_interfaces
      implicit none

      type(user_input_parameters)  :: param_dump !! Data type of dumped parameter file
      integer(I4B), parameter :: LUN = 7 !! Unit number of output file
      integer(I4B)            :: ierr     !! Error code
      integer(I4B), save      :: idx = 1  !! Index of current dump file. Output flips between 2 files for extra security
                                          !!    in case the program halts during writing
      character(STRMAX)       :: error_message           !! Error message in UDIO procedure

      param_dump = param
      param_dump%t0 = t
      param_dump%inplfile = trim(adjustl(DUMP_PL_FILE(idx))) 
      param_dump%intpfile = trim(adjustl(DUMP_TP_FILE(idx))) 
      open(unit = LUN, file = DUMP_PARAM_FILE(idx), status='replace', form = 'formatted', iostat =ierr)
      if (ierr /=0) then
         write(*,*) 'Swiftest error.'
         write(*,*) '   Could not open dump file: ',trim(adjustl(DUMP_PARAM_FILE(idx)))
         call util_exit(FAILURE)
      end if
      
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') param_dump
      call param_dump%udio_writer(LUN, iotype="none",v_list=(/0/),iostat=ierr,iomsg=error_message)

      if (idx == 2) then
         idx = 1
      else
         idx = 2
      end if

      close(LUN)

      return

      end procedure user_dump_param
end submodule s_user_dump_param
