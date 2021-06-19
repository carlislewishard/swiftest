submodule(user) s_user_udio_writer
contains
   module subroutine user_udio_writer(param, unit, iotype, v_list, iostat, iomsg) 
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      use swiftest, except_this_one => user_udio_writer
      implicit none

      ! Arguments
      class(user_input_parameters),intent(in)  :: param         !! Output collection of user-defined parameters
      integer, intent(in)                 :: unit        
      character(len=*), intent(in)        :: iotype
      integer, intent(in)                 :: v_list(:)
      integer, intent(out)                :: iostat
      character(len=*), intent(inout)     :: iomsg

      ! Internals
      character(*),parameter :: Ifmt  = '(I0)'         !! Format label for integer values
      character(*),parameter :: Rfmt  = '(ES25.17)'    !! Format label for real values 
      character(*),parameter :: Rarrfmt  = '(3(ES25.17,1X))'    !! Format label for real values 
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 
      character(len=*), parameter :: Afmt = '(A25,1X,64(:,A25,1X))'
      character(256)          :: param_name, param_value
      type character_array
         character(25) :: value
      end type character_array
      type(character_array), dimension(:), allocatable :: param_array
      integer(I4B) :: i

      write(param_name, Afmt) "T0"; write(param_value,Rfmt) param%t0; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "TSTOP"; write(param_value, Rfmt) param%tstop; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "DT"; write(param_value, Rfmt) param%dt; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "PL_IN"; write(param_value, Afmt) trim(adjustl(param%inplfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "TP_in"; write(param_value, Afmt) trim(adjustl(param%intpfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "IN_TYPE"; write(param_value, Afmt) trim(adjustl(param%in_type)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      if (param%istep_out > 0) then
         write(param_name, Afmt) "ISTEP_OUT"; write(param_value, Ifmt) param%istep_out; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "BIN_OUT"; write(param_value, Afmt) trim(adjustl(param%outfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "PARTICLE_FILE"; write(param_value, Afmt) trim(adjustl(param%particle_file)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "OUT_TYPE"; write(param_value, Afmt) trim(adjustl(param%out_type)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "OUT_FORM"; write(param_value, Afmt) trim(adjustl(param%out_form)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "OUT_STAT"; write(param_value, Afmt) "APPEND"; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      end if
      write(param_name, Afmt) "ENC_OUT"; write(param_value, Afmt) trim(adjustl(param%encounter_file)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      if (param%istep_dump > 0) then
         write(param_name, Afmt) "ISTEP_DUMP"; write(param_value, Ifmt) param%istep_dump; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      end if
      if (param%j2rp2 > VSMALL) then
         write(param_name, Afmt) "J2 "; write(param_value, Rfmt) param%j2rp2; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         if (param%j4rp4 > VSMALL) then
            write(param_name, Afmt) "J4 "; write(param_value, Rfmt) param%j4rp4; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         end if
      end if
      write(param_name, Afmt) "CHK_RMIN"; write(param_value, Rfmt) param%rmin; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "CHK_RMAX"; write(param_value, Rfmt) param%rmax; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "CHK_EJECT"; write(param_value, Rfmt) param%rmaxu; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "CHK_QMIN"; write(param_value, Rfmt) param%qmin; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      if (param%qmin >= 0.0_DP) then
         write(param_name, Afmt) "CHK_QMIN_COORD"; write(param_value, Afmt) trim(adjustl(param%qmin_coord)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         allocate(param_array(2))
         write(param_array(1)%value, Rfmt) param%qmin_alo
         write(param_array(2)%value, Rfmt) param%qmin_ahi
         write(param_name, Afmt) "CHK_QMIN_RANGE"; write(unit, Afmt) adjustl(param_name), adjustl(param_array(1)%value), adjustl(param_array(2)%value)
      end if
      if (param%lmtiny) then
         write(param_name, Afmt) "MTINY"; write(param_value, Rfmt) param%mtiny; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      end if
      write(param_name, Afmt) "MU2KG"; write(param_value, Rfmt) MU2KG; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "TU2S"; write(param_value, Rfmt) TU2S ; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "DU2M"; write(param_value, Rfmt) DU2M; write(unit, Afmt) adjustl(param_name), adjustl(param_value)

      ! Special handling is required for writing the random number seed array as its size is not known until runtime
      ! For the "SEED" parameter line, the first value will be the size of the seed array and the rest will be the seed array elements
      write(param_name, Afmt) "SEED"
      if (allocated(param_array)) deallocate(param_array)
      allocate(param_array(0:size(param%seed)))
      write(param_array(0)%value, Ifmt) size(param%seed)
      do i = 1, size(param%seed)
         write(param_array(i)%value, Ifmt) param%seed(i)
      end do
      write(unit, Afmt, advance='no') adjustl(param_name), adjustl(param_array(0)%value)
      do i = 1, size(param%seed)
         if (i < size(param%seed)) then
            write(unit, Afmt, advance='no') adjustl(param_array(i)%value)
         else
            write(unit, Afmt) adjustl(param_array(i)%value)
         end if
      end do
      ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      write(param_name, Afmt) "EXTRA_FORCE"; write(param_value, Lfmt) param%lextra_force; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "BIG_DISCARD"; write(param_value, Lfmt) param%lbig_discard; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "RHILL_PRESENT"; write(param_value, Lfmt) param%lrhill_present; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "CHK_CLOSE"; write(param_value, Lfmt) param%lclose; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "FRAGMENTATION"; write(param_value, Lfmt)  param%lfragmentation; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "ENERGY"; write(param_value, Lfmt)  param%lenergy; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "ROTATION"; write(param_value, Lfmt)  param%lrotation; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "TIDES"; write(param_value, Lfmt)  param%ltides; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "GR"; write(param_value, Lfmt)  param%lgr; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "YARKOVSKY"; write(param_value, Lfmt)  param%lyarkovsky; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "YORP"; write(param_value, Lfmt)  param%lyorp; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      write(param_name, Afmt) "RINGMOONS"; write(param_value, Lfmt)  param%lringmoons; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      if (param%lringmoons) then
         write(param_name, Afmt) "RING_OUTFILE"; write(param_value, Afmt) trim(adjustl(param%ring_outfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
      end if
      if (param%lenergy) then
         write(param_name, Afmt) "FIRSTENERGY"; write(param_value, Lfmt) param%lfirstenergy; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "EORBIT_ORIG"; write(param_value, Rfmt) param%Eorbit_orig; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "MTOT_ORIG"; write(param_value, Rfmt) param%Mtot_orig; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(unit, '("LTOT_ORIG  ",3(1X,ES25.17))') param%Ltot_orig(:)
         write(unit, '("LORBIT_ORIG",3(1X,ES25.17))') param%Lorbit_orig(:)
         write(unit, '("LSPIN_ORIG ",3(1X,ES25.17))') param%Lspin_orig(:)
      end if
      write(param_name, Afmt) "FIRSTKICK"; write(param_value, Lfmt) param%lfirstkick; write(unit, Afmt) adjustl(param_name), adjustl(param_value)

      iostat = 0

      return

   end subroutine user_udio_writer
end submodule s_user_udio_writer
