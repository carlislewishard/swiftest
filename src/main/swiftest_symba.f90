program swiftest_symba
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for the Symplectic Massive Body Algorithm
   !!
   !! Adapted from Swifter by David E. Kaufmanna swiftert_symba.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift program swift_symba5.f
   !! Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
   use swiftest

   !> The following are temporary until the conversion to the new module structure is complete
   use module_symba
   use module_interfaces
   use module_swiftestalloc
   implicit none

   ! Arguments
   type(user_input_parameters)  :: param    ! derived type containing user-defined parameters
   integer(I4B)      :: istep_out      ! time steps between binary outputs
   integer(I4B)      :: istep_dump     ! time steps between dumps
   real(DP)          :: t0             ! integration start time
   real(DP)          :: tstop          ! integration stop time
   real(DP)          :: dt             ! time step
   real(DP)          :: j2rp2          ! j2*r^2 term for central body
   real(DP)          :: j4rp4          ! j4*r^4 term for central body
   real(DP)          :: rmin           ! minimum heliocentric radius for test particle
   real(DP)          :: rmax           ! maximum heliocentric radius for test particle
   real(DP)          :: rmaxu          ! maximum unbound heliocentric radius for test particle
   real(DP)          :: qmin           ! minimum pericenter distance for test particle
   real(DP)          :: qmin_alo       ! minimum semimajor axis for qmin
   real(DP)          :: qmin_ahi       ! maximum semimajor axis for qmin
   character(strmax) :: qmin_coord     ! coordinate frame to use for qmin
   character(strmax) :: encounter_file ! name of output file for encounters
   character(strmax) :: inplfile       ! name of input file for planets
   character(strmax) :: intpfile       ! name of input file for test particles
   character(strmax) :: in_type        ! format of input data files
   character(strmax) :: outfile        ! name of output binary file
   character(strmax) :: out_type       ! binary format of output file
   character(strmax) :: out_form       ! data to write to output file
   character(strmax) :: out_stat       ! open status for output binary file

   ! Internals
   logical                       :: lfrag_add, ldiscard, ldiscard_tp
   integer(I4B)                  :: npl, nplm, ntp, ntp0, nsppl, nsptp, iout, idump, iloop
   integer(I4B)                  :: nplplenc, npltpenc, nmergeadd, nmergesub
   real(DP)                      :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset, Loffset, msys
   real(DP), dimension(NDIM)     :: htot
   real(DP)                      :: te_orig, te_error, te_off_error, te_after, te_before
   real(DP)                      :: Mtot_orig, Mtot_now, Merror
   real(DP)                      :: Ltot_orig, Ltot_now, Lerror, L_off_error
   character(STRMAX)             :: inparfile
   type(symba_pl)                :: symba_plA
   type(symba_tp)                :: symba_tpA
   type(swiftest_tp)             :: discard_tpA
   type(swiftest_pl)             :: discard_plA
   type(symba_plplenc)           :: plplenc_list
   type(symba_pltpenc)           :: pltpenc_list
   type(symba_merger)            :: mergeadd_list, mergesub_list
   integer(I4B), parameter       :: egyiu = 72
   real(DP)                      :: start, finish
   INTEGER(I8B)                  :: clock_count, count_rate, count_max
   INTEGER(I4B)                  :: ierr
   INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: k_plpl, k_pltp
   INTEGER(I8B)                  :: num_plpl_comparisons, num_pltp_comparisons
   

! Executable code

   call util_version

   call get_command_argument(1, inparfile, status = ierr) 
   if (ierr /= 0) then
       write(*, 100, advance = "NO") "Enter name of parameter data file: "
       read(*, 100) inparfile
   end if
   
   100 format(a)
   inparfile = trim(adjustl(inparfile))
   ! read in the param.in file and get simulation parameters
   call param%read_from_file(inparfile)
   param%lmtiny = .true. ! Turn this on for SyMBA

   ! temporary until the conversion to the derived type argument list is complete
   t0 = param%t0
   tstop = param%tstop
   dt = param%dt
   inplfile = param%inplfile
   intpfile = param%intpfile
   in_type = param%in_type
   istep_out = param%istep_out
   outfile = param%outfile
   out_type = param%out_type
   out_form = param%out_form
   out_stat = param%out_stat
   istep_dump = param%istep_dump
   j2rp2 = param%j2rp2
   j4rp4 = param%j4rp4
   rmin = param%rmin
   rmax = param%rmax
   rmaxu = param%rmaxu
   qmin = param%qmin
   qmin_coord = param%qmin_coord
   qmin_alo = param%qmin_alo
   qmin_ahi = param%qmin_ahi
   encounter_file = param%encounter_file
   mtiny = param%mtiny
   !^^^^^^^^^^^^^^^^^^^^^^^^^
   if (.not. param%lrhill_present) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Integrator SyMBA requires massive body Hill sphere radii on input"
      call util_exit(failure)
   end if

   ! reads in initial conditions of all massive bodies from input file
   call symba_plA%helio%swiftest%read_from_file(param)
   call symba_tpA%helio%swiftest%read_from_file(param)

   !Temporary until the argument lists get fixed
      npl = symba_plA%helio%swiftest%nbody
      ntp = symba_tpA%helio%swiftest%nbody
   ! Temporary fix until all of the data structures are converted to OOP and inheritance works properly
      call symba_plA%helio%swiftest%dealloc()
      call symba_tpA%helio%swiftest%dealloc()
      call symba_pl_allocate(symba_plA,npl)
      call symba_tp_allocate(symba_tpA,ntp)
      call symba_plA%helio%swiftest%dealloc()
      call symba_tpA%helio%swiftest%dealloc()
      call symba_plA%helio%swiftest%read_from_file(param)
      call symba_tpA%helio%swiftest%read_from_file(param)
   !**************************************************

   if (ntp > 0) then
      call symba_pltpenc_allocate(pltpenc_list, ntp)
   end if

   ! reorder by mass 
   call symba_reorder_pl(npl, symba_plA)
   call util_valid(npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest)
   ntp0 = ntp
   t = t0
   tbase = t0
   iloop = 0
   iout = istep_out
   idump = istep_dump
   nmergeadd = 0
   nmergesub = 0
   nsppl = 0
   nsptp = 0
   if ((istep_out > 0).and.(out_stat == "NEW")) then
      call io_write_frame(t, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, param) 
   end if
   if (out_stat == "OLD") then
      open(unit = egyiu, file = energy_file, form = "formatted", status = "old", action = "write", position = "append")
   else 
      open(unit = egyiu, file = energy_file, form = "formatted", status = "replace", action = "write")
   end if
   300 format(10(1x, e23.16))
   nplm = count(symba_plA%helio%swiftest%mass>mtiny)
   CALL util_dist_index_plpl(npl, nplm, num_plpl_comparisons, k_plpl)
   CALL util_dist_index_pltp(nplm, ntp, num_pltp_comparisons, k_pltp)

   if (param%lenergy) then
      eoffset = 0.0_DP
      Loffset = 0.0_DP
      if(num_plpl_comparisons > param%eucl_threshold) then
         call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, te, htot, msys)
      else
         call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
      end if
      Ltot_orig = NORM2(htot)
      te_orig = te
      Mtot_orig = msys
      write(egyiu,300) t, ke, pe, te, htot, eoffset, msys
   end if
   write(*, *) " *************** Main Loop *************** "

   call system_clock(clock_count, count_rate, count_max)
   start = clock_count / (count_rate * 1.0_DP)
   do while ((t < tstop) .and. ((ntp0 == 0) .or. (ntp > 0)))
      if(num_plpl_comparisons > param%eucl_threshold) then
         call symba_step_eucl(t, dt, param,npl,ntp,symba_plA, symba_tpa,nplplenc, npltpenc,&
               plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, Loffset, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
      else
         call symba_step(t, dt, param,npl,ntp,symba_plA, symba_tpA, nplplenc, npltpenc,&
               plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, Loffset)
      end if
      ldiscard = .false. 
      ldiscard_tp = .false.
      lfrag_add = .false.

      if (param%lenergy) then
         if(num_plpl_comparisons > param%eucl_threshold) then
            call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, &
                  ke, pe, te, htot, msys)
         else
            call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
         end if
         te_before = te
      end if

      call symba_discard_merge_pl(symba_plA, nplplenc, plplenc_list, ldiscard, mergeadd_list, nmergeadd)                                  
      call symba_discard_pl(t, npl, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, ldiscard)
      call symba_discard_tp(t, npl, ntp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &    
            qmin_alo, qmin_ahi, param%lrhill_present, ldiscard_tp)
      if (ldiscard .or. ldiscard_tp .or. lfrag_add) then
         call symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
            discard_tpA, param, ldiscard, ldiscard_tp)

         if (ldiscard .or. ldiscard_tp) then
            call io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergesub, symba_plA, &
               discard_plA, discard_tpA, mergeadd_list, mergesub_list, discard_file, param%lbig_discard) 
            nmergeadd = 0
            nmergesub = 0
            nsppl = 0
            nsptp = 0
            deallocate(k_plpl)
            nplm = count(symba_plA%helio%swiftest%mass(1:npl) > mtiny)
            CALL util_dist_index_plpl(npl, nplm, num_plpl_comparisons, k_plpl)
            if (ntp > 0) then
                 deallocate(k_pltp)
                 call util_dist_index_pltp(nplm, ntp, num_pltp_comparisons, k_pltp)          
            end if 
         end if

         if (param%lenergy) then
            if(num_plpl_comparisons > param%eucl_threshold) then
               call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, &
                  ke, pe, te, htot, msys)
            else
               call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
            end if
            te_after = te
            te_error = (te_after - te_orig) / te_orig     ! Total energy error of system now compared to original system
            eoffset = eoffset + (te_after - te_before)    ! Total running energy offset from collision in this step
            te_off_error = te_error + (eoffset / te_orig) ! Total energy of the system plus any energy lost due to collisions
            write(egyiu,300) t, ke, pe, te, htot, eoffset, Loffset, msys
         end if
      end if


      !if (param%lenergy) then
      !   if(num_plpl_comparisons > param%eucl_threshold) then
      !      call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, &
      !            ke, pe, te, htot, msys)
      !   else
      !      call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
      !   end if
      !   te_after = te
      !   te_error = (te_after - te_orig) / te_orig     ! Total energy error of system now compared to original system
      !   eoffset = eoffset + (te_after - te_before)    ! Total running energy offset from collision in this step
      !   te_off_error = te_error + (eoffset / te_orig) ! Total energy of the system plus any energy lost due to collisions

         !write(*, *) "  DE/E0 = ", te_error, "; eoffset = ", eoffset, "; (DE+eoffset)/E0 = ", te_off_error
      !end if

      iloop = iloop + 1
      if (iloop == loopmax) then
          tbase = tbase + iloop*dt
          iloop = 0
      end if
      t = tbase + iloop*dt

      if (istep_out > 0) then
         iout = iout - 1
         if (iout == 0) then
            call io_write_frame(t, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, param)
            iout = istep_out
            if (param%lenergy) then
               if(num_plpl_comparisons > param%eucl_threshold) then
                  call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, &
                     te, htot, msys)
               else
                  call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
               end if
               write(egyiu,300) t, ke, pe, te, htot, eoffset, Loffset, msys
            end if
         end if
      end if

      if (istep_dump > 0) then
         idump = idump - 1
         if (idump == 0) then
            if (param%lenergy) then
               if(num_plpl_comparisons > param%eucl_threshold) then
                  call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, &
                     te, htot, msys)
               else
                  call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
               end if
               Mtot_now = msys
               Merror = (Mtot_now - Mtot_orig) / Mtot_orig
               Ltot_now = norm2(htot)
               Lerror = (Ltot_now - Ltot_orig) / Ltot_orig
               L_off_error = (Ltot_now - Ltot_orig + Loffset) / Ltot_orig
               te_error = (te - te_orig) / te_orig
               te_off_error = (te + eoffset - te_orig) / te_orig
            end if
            tfrac = (t - t0)/(tstop - t0)
            write(*, 200) t, tfrac, npl, ntp
200         format(" Time = ", es12.5, "; fraction done = ", f5.3, "; number of active pl, tp = ", i7, ", ", i7)

            call system_clock(clock_count)
            finish = clock_count / (count_rate * 1.0_DP)
            write(*,*) "      Wall time (s): ", finish - start

205         format("  DL/L0 = ", ES12.5, "; (DL+Loffset)/L0 = ", ES12.5, &
                   "; DE/E0 = ", ES12.5, "; (DE+eoffset)/E0 = ", ES12.5, &
                   "; DM/M0 = ", ES12.5)
            if (param%lenergy) write(*, 205) Lerror, L_off_error, te_error, te_off_error, Merror
            call param%dump_to_file(t)
            call io_dump_pl(npl, symba_plA%helio%swiftest, param%lclose, param%lrhill_present)
            call io_dump_tp(ntp, symba_tpA%helio%swiftest)
            idump = istep_dump
         end if
      end if

      if (allocated(plplenc_list%status)) then
         plplenc_list%status(:) = 0
         plplenc_list%lvdotr(:) = .false.
         plplenc_list%level(:) = 0
         plplenc_list%index1(:) = 0
         plplenc_list%index2(:) = 0
         plplenc_list%enc_child(:) = 0 
         plplenc_list%enc_parent(:) = 0
      end if

      if (allocated(pltpenc_list%status)) then
         pltpenc_list%status(:) = 0
         pltpenc_list%lvdotr(:) = .false.
         pltpenc_list%level(:) = 0
         pltpenc_list%indexpl(:) = 0
         pltpenc_list%indextp(:) = 0
      end if

      if (allocated(mergeadd_list%name)) then
         mergeadd_list%name(:) = 0
         mergeadd_list%index_ps(:) = 0
         mergeadd_list%status(:) = 0
         mergeadd_list%ncomp(:) = 0
         mergeadd_list%xh(:,:) = 0
         mergeadd_list%vh(:,:) = 0
         mergeadd_list%mass(:) = 0
         mergeadd_list%radius(:) = 0
         mergeadd_list%IP(:,:) = 0
         mergeadd_list%rot(:,:) = 0
      end if

      if (allocated(mergesub_list%name)) then
         mergesub_list%name(:) = 0
         mergesub_list%index_ps(:) = 0
         mergesub_list%status(:) = 0
         mergesub_list%ncomp(:) = 0
         mergesub_list%xh(:,:) = 0
         mergesub_list%vh(:,:) = 0
         mergesub_list%mass(:) = 0
         mergesub_list%radius(:) = 0
         mergeadd_list%IP(:,:) = 0
         mergeadd_list%rot(:,:) = 0
      end if

      if (allocated(discard_plA%name)) call swiftest_pl_deallocate(discard_plA)
      if (allocated(discard_tpA%name)) call swiftest_tp_deallocate(discard_tpA)

   end do

   if (param%lenergy) then
      if(num_plpl_comparisons > param%eucl_threshold) then
         call symba_energy_eucl(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, k_plpl, num_plpl_comparisons, ke, pe, te, htot, msys)
      else
         call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot, msys)
      end if
      Mtot_now = msys
      Merror = (Mtot_now - Mtot_orig) / Mtot_orig
      Ltot_now = norm2(htot)
      Lerror = (Ltot_now - Ltot_orig) / Ltot_orig
      L_off_error = (Ltot_now - Ltot_orig + Loffset) / Ltot_orig
      te_error = (te - te_orig) / te_orig
      te_off_error = (te + eoffset - te_orig) / te_orig
      write(*,*) 'Final angular momentum and energy errors'
      write(*, 205) Lerror, L_off_error, te_error, te_off_error, Merror
      write(egyiu,300) t, ke, pe, te, htot, eoffset, Loffset, msys
      close(egyiu)
   end if

   call param%dump_to_file(t)
   call io_dump_pl(npl, symba_plA%helio%swiftest, param%lclose, param%lrhill_present)
   call io_dump_tp(ntp, symba_tpA%helio%swiftest)

   call symba_pl_deallocate(symba_plA)
   call symba_merger_deallocate(mergeadd_list)
   call symba_merger_deallocate(mergesub_list)
   call symba_plplenc_deallocate(plplenc_list)
   if (ntp > 0) then
      call symba_tp_deallocate(symba_tpA)
      call symba_pltpenc_deallocate(pltpenc_list)
   end if
   call system_clock(clock_count)
   finish = clock_count / (count_rate * 1.0_DP)
   write(*,*) 'Wall time to complete run (s): ', finish - start

   call util_exit(SUCCESS)

   stop

end program swiftest_symba

