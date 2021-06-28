program swiftest_symba
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for the Symplectic Massive Body Algorithm
   !!
   !! Adapted from Swifter by David E. Kaufmanna swifter_symba.f90
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

   ! Internals
   logical                       :: lfrag_add, ldiscard_pl, ldiscard_tp
   integer(I4B)                  :: nplm, ntp, ntp0, nsppl, nsptp, iout, idump, iloop, i
   integer(I4B)                  :: nplplenc, npltpenc, nmergeadd, nmergesub
   real(DP)                      :: t, tfrac, tbase,  msys
   real(DP)                      :: Ecollision, Eorbit_before, Eorbit_after, ke_orbit_before, ke_spin_before, ke_orbit_after, ke_spin_after, pe_before, pe_after
   real(DP), dimension(NDIM)     :: Lorbit, Lspin
   character(STRMAX)             :: inparfile
   type(symba_pl)                :: symba_plA
   type(symba_tp)                :: symba_tpA
   type(symba_tp)                :: discard_tpA
   type(symba_pl)                :: discard_plA
   type(symba_plplenc)           :: plplenc_list
   type(symba_pltpenc)           :: pltpenc_list
   type(symba_merger)            :: mergeadd_list, mergesub_list
   real(DP)                      :: start, finish, deltawall, wallperstep
   integer(I8B)                  :: clock_count, count_rate, count_max
   integer(I4B)                  :: ierr
   integer(I4B)                  :: nplfile, ntpfile ! Temporary variables until we switch over to the OOP branch
   integer(I4B), parameter       :: EGYDUMP = 88
   character(len=*), parameter   :: simtimefmt = '(" Time = ", es12.5, "; fraction done = ", f5.3, "; number of active pl, tp = ", i7, ", ", i7)'
   character(len=*), parameter   :: walltimefmt = '("      Wall time (s): ", es12.5, "; Wall time/step in this interval (s):  ", es12.5)'
   character(len=*), parameter   :: endwallfmt = '("Wall time to complete run (s): ", es12.5)'
   integer(I4B), dimension(:), allocatable :: discard_stat_list
   logical, dimension(:), allocatable :: discard_l_pl

! Executable code
   ! temporary until the conversion to the derived type argument list is complete
   associate(  t0 => param%t0, &
               tstop => param%tstop, &
               dt => param%dt, &
               inplfile => param%inplfile, &
               intpfile => param%intpfile, &
               in_type => param%in_type, &
               istep_out => param%istep_out, &
               outfile => param%outfile, &
               out_type => param%out_type, &
               out_form => param%out_form, &
               out_stat => param%out_stat, &
               istep_dump => param%istep_dump, &
               j2rp2 => param%j2rp2, &
               j4rp4 => param%j4rp4, &
               rmin => param%rmin, &
               rmax => param%rmax, &
               rmaxu => param%rmaxu, &
               qmin => param%qmin, &
               qmin_coord => param%qmin_coord, &
               qmin_alo => param%qmin_alo, &
               qmin_ahi => param%qmin_ahi, &
               encounter_file => param%encounter_file, &
               mtiny => param%mtiny, &
               npl => symba_plA%helio%swiftest%nbody, &
               ntp => symba_tpA%helio%swiftest%nbody)

      call util_version
       
      call get_command_argument(1, inparfile, status = ierr) 
      if (ierr /= 0) then
         write(*, 100, advance = "NO") "Enter name of parameter data file: "
         read(*, 100) inparfile
      end if
      
      100 format(a)
      inparfile = trim(adjustl(inparfile))
      ! read in the param.in file and get simulation parameters
      call param%read_from_file(inparfile, symba_plA%helio%swiftest)
      param%lmtiny = .true. ! Turn this on for SyMBA
      
      ! reads in initial conditions of all massive bodies from input file
      call symba_read_pl_in(symba_plA, param) 
      call symba_tpA%helio%swiftest%read_from_file(param)

      !Temporary until the argument lists get fixed
      ! Temporary fix until all of the data structures are converted to OOP and inheritance works properly
         nplfile = npl
         ntpfile = ntp
         call symba_plA%helio%swiftest%dealloc()
         call symba_tpA%helio%swiftest%dealloc()
         call symba_pl_allocate(symba_plA,nplfile)
         call symba_tp_allocate(symba_tpA,ntpfile)
         call symba_plA%helio%swiftest%dealloc()
         call symba_tpA%helio%swiftest%dealloc()
         call symba_read_pl_in(symba_plA, param)
         call symba_tpA%helio%swiftest%read_from_file(param)
      !**************************************************
      call io_write_particle_pl(symba_plA%helio%swiftest, [(i, i=1,npl)], param)

      if (ntp > 0) then
         call symba_pltpenc_allocate(pltpenc_list, ntp)
      end if

      ! reorder by mass 
      if (out_stat /= "APPEND") call symba_reorder_pl(npl, symba_plA)  ! This is a new run, so we will sort the massive body list by mass
          
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
      if ((istep_out > 0).and.((out_stat == "NEW").or.(out_stat == "REPLACE"))) then
         call io_write_frame(t, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, param) 
      end if

      nplm = count(symba_plA%helio%swiftest%mass(1:npl) > mtiny)
      CALL util_dist_index_plpl(npl, nplm, symba_plA)
      !CALL util_dist_index_pltp(nplm, ntp, symba_tpA)

      ! Save initial mass and angular momentum of the central body
      symba_plA%helio%swiftest%Mcb_initial = symba_plA%helio%swiftest%mass(1)
      symba_plA%helio%swiftest%Rcb_initial = symba_plA%helio%swiftest%radius(1)
      symba_plA%helio%swiftest%Lcb_initial(:) = symba_plA%helio%swiftest%Ip(3,1) * symba_plA%helio%swiftest%mass(1) * &
                                                symba_plA%helio%swiftest%radius(1)**2 * symba_plA%helio%swiftest%rot(:,1)

      if (param%lenergy) then
         call coord_h2b(npl, symba_plA%helio%swiftest, msys)
         call io_conservation_report(t, symba_plA, npl, j2rp2, j4rp4, param, lterminal=.false.) 
         call symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke_orbit_before, ke_spin_before, pe_before, Lorbit, Lspin)

      end if
      write(*, *) " *************** Main Loop *************** "

      call system_clock(clock_count, count_rate, count_max)
      start = clock_count / (count_rate * 1.0_DP)
      finish = start
      do while ((t < tstop) .and. ((ntp0 == 0) .or. (ntp > 0)))
         call util_hills(npl, symba_plA%helio%swiftest)
         call symba_step_eucl(t, dt, param,npl,ntp,symba_plA, symba_tpA, nplplenc, npltpenc,&
               plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list)

         if (allocated(discard_l_pl)) deallocate(discard_l_pl)
         allocate(discard_l_pl(npl))
         discard_l_pl(:) = .false.
         ldiscard_tp = .false.
         lfrag_add = .false.
         call symba_discard_pl(t, npl, ntp, symba_plA, symba_tpA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, discard_l_pl, discard_stat_list)
         ldiscard_pl = any(discard_l_pl(:))
         call symba_discard_tp(t, npl, ntp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &    
                                qmin_alo, qmin_ahi, param%lrhill_present, ldiscard_tp)
         call symba_collision(t, symba_plA, nplplenc, plplenc_list, lfrag_add, mergeadd_list, nmergeadd, param)
         ldiscard_pl = ldiscard_pl .or. lfrag_add

         if (ldiscard_pl .or. ldiscard_tp) then
            if (param%lenergy) then
               call symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke_orbit_before, ke_spin_before, pe_before, Lorbit, Lspin)
               Eorbit_before = ke_orbit_before + ke_spin_before + pe_before
            end if
            call symba_rearray(t, npl, nplm, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
                               discard_tpA, ldiscard_pl, ldiscard_tp, mtiny, param, discard_l_pl, discard_stat_list)
            call io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergesub, symba_plA, &
                                        discard_plA%helio%swiftest, discard_tpA%helio%swiftest, mergeadd_list, mergesub_list, discard_file, param%lbig_discard) 
            nmergeadd = 0
            nmergesub = 0
            nsppl = 0
            nsptp = 0
            nplm = count(symba_plA%helio%swiftest%mass(1:npl) > mtiny)

            if (param%lenergy)  then
               call symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke_orbit_after, ke_spin_after, pe_after, Lorbit, Lspin)
               Eorbit_after = ke_orbit_after + ke_spin_after + pe_after
               Ecollision = Eorbit_after - Eorbit_before   ! Energy change resulting in this collisional event Total running energy offset from collision in this step
               symba_plA%helio%swiftest%Ecollisions = symba_plA%helio%swiftest%Ecollisions + Ecollision
            end if
            !if (ntp > 0) call util_dist_index_pltp(nplm, ntp, symba_plA, symba_tpA)
         end if

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
            end if
         end if

         if (istep_dump > 0) then
            idump = idump - 1
            if (idump == 0) then
               if (param%lenergy) call io_conservation_report(t, symba_plA, npl, j2rp2, j4rp4, param, lterminal=.true.) 
               tfrac = (t - t0) / (tstop - t0)
               write(*, simtimefmt) t, tfrac, npl, ntp

               call system_clock(clock_count)
               deltawall = clock_count / (count_rate * 1.0_DP) - finish
               wallperstep = deltawall / istep_dump
               finish = clock_count / (count_rate * 1.0_DP)
               write(*,walltimefmt) finish - start, wallperstep

               call param%dump_to_file(t, symba_plA%helio%swiftest)
               call io_dump_pl_symba(npl, symba_plA, param)
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
            plplenc_list%xh1(:,:) = 0
            plplenc_list%xh2(:,:) = 0
            plplenc_list%vb1(:,:) = 0
            plplenc_list%vb2(:,:) = 0
         end if

         if (allocated(pltpenc_list%status)) then
            pltpenc_list%status(:) = 0
            pltpenc_list%lvdotr(:) = .false.
            pltpenc_list%level(:) = 0
            pltpenc_list%indexpl(:) = 0
            pltpenc_list%indextp(:) = 0
         end if

         if (allocated(mergeadd_list%id)) then
            mergeadd_list%id(:) = 0
            mergeadd_list%index_ps(:) = 0
            mergeadd_list%status(:) = 0
            mergeadd_list%xb(:,:) = 0
            mergeadd_list%vb(:,:) = 0
            mergeadd_list%mass(:) = 0
            mergeadd_list%radius(:) = 0
            mergeadd_list%Ip(:,:) = 0
            mergeadd_list%rot(:,:) = 0
         end if

         if (allocated(mergesub_list%id)) then
            mergesub_list%id(:) = 0
            mergesub_list%index_ps(:) = 0
            mergesub_list%status(:) = 0
            mergesub_list%xb(:,:) = 0
            mergesub_list%vb(:,:) = 0
            mergesub_list%mass(:) = 0
            mergesub_list%radius(:) = 0
            mergeadd_list%Ip(:,:) = 0
            mergeadd_list%rot(:,:) = 0
         end if

         if (allocated(discard_plA%helio%swiftest%id)) call symba_pl_deallocate(discard_plA)
         if (allocated(discard_tpA%helio%swiftest%id)) call symba_tp_deallocate(discard_tpA)

      end do

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
      write(*, endwallfmt) finish - start

      call util_exit(SUCCESS)
   end associate

   stop

end program swiftest_symba

