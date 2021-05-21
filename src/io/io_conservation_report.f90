submodule (io) s_io_conservation_report
contains
   module subroutine io_conservation_report(t, symba_plA, npl, j2rp2, j4rp4, param, lterminal)
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Reports the current state of energy, mass, and angular momentum conservation in a run
      use module_interfaces
      use module_symba
      implicit none
      ! Arguments
      real(DP),                    intent(in)    :: t            !! Current time of simulation
      type(symba_pl),              intent(inout) :: symba_plA    !! Swiftest planet data structure
      integer(I4B),                intent(in)    :: npl          !! Number of massive bodies
      real(DP),                    intent(in)    :: j2rp2, j4rp4 !! Central body oblateness terms
      type(user_input_parameters), intent(in)    :: param        !! Input colleciton of user-defined parameters
      logical,                     intent(in)    :: lterminal    !! Indicates whether to output information to the terminal screen
      ! Internals
      real(DP), dimension(NDIM), save :: Ltot_orig, Ltot_last
      real(DP), save                  :: Eorbit_orig, Mtot_orig, Lmag_orig, ke_orb_last, ke_spin_last, pe_last, Eorbit_last
      real(DP)                        :: ke_orbit, ke_spin, pe, Eorbit
      real(DP), dimension(NDIM)       :: Ltot_now
      real(DP)                        :: Eorbit_error, Etotal_error, Ecoll_error
      real(DP)                        :: Mtot_now, Merror
      real(DP)                        :: Lmag_now, Lerror
      logical, save                   :: lfirst = .true.
      character(len=*), parameter     :: egyfmt = '(ES23.16,10(",",ES23.16,:))' ! Format code for all simulation output
      character(len=*), parameter     :: egyheader = '("t,Eorbit,Ecollisions,Lx,Ly,Lz,Mtot")'
      integer(I4B), parameter         :: egyiu = 72
      character(len=*), parameter     :: egytermfmt = '("  DL/L0 = ", ES12.5 &
                                                         "; DEcollisions/|E0| = ", ES12.5, &
                                                         "; D(Eorbit+Ecollisions)/|E0| = ", ES12.5, &
                                                         "; DM/M0 = ", ES12.5)'

      associate(Ecollisions => symba_plA%helio%swiftest%Ecollisions, Lescape => symba_plA%helio%swiftest%Lescape, Mescape => symba_plA%helio%swiftest%Mescape, &
                  mass => symba_plA%helio%swiftest%mass, dMcb => symba_plA%helio%swiftest%dMcb, Mcb_initial => symba_plA%helio%swiftest%Mcb_initial)
         if (lfirst) then
            if (param%out_stat == "OLD") then
               open(unit = egyiu, file = energy_file, form = "formatted", status = "old", action = "write", position = "append")
            else 
               open(unit = egyiu, file = energy_file, form = "formatted", status = "replace", action = "write")
               write(egyiu,egyheader)
            end if
         end if
         call symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke_orbit, ke_spin, pe, Eorbit, Ltot_now)
         Mtot_now = dMcb + sum(mass(2:npl)) + Mcb_initial + Mescape
         Ltot_now(:) = Lescape(:) + Ltot_now(:)
         if (lfirst) then
            Eorbit_orig = Eorbit
            Mtot_orig = Mtot_now
            Ltot_orig(:) = Ltot_now(:)
            Lmag_orig = norm2(Ltot_orig(:))
            lfirst = .false.
         end if

         write(egyiu,egyfmt) t, Eorbit, Ecollisions, Ltot_now, Mtot_now
         flush(egyiu)
         if (.not.lfirst .and. lterminal) then 
            Lmag_now = norm2(Ltot_now)
            Lerror = norm2(Ltot_now - Ltot_orig) / Lmag_orig
            Eorbit_error = (Eorbit - Eorbit_orig) / abs(Eorbit_orig)
            Ecoll_error = -Ecollisions / abs(Eorbit_orig)
            Etotal_error = (Eorbit - (Eorbit_orig - Ecollisions)) / abs(Eorbit_orig)
            Merror = (Mtot_now - Mtot_orig) / Mtot_orig
            write(*, egytermfmt) Lerror, Ecoll_error, Etotal_error, Merror
            if (Ecoll_error > 0.0_DP) then
               write(*,*) 'Something has gone wrong! Collisional energy should not be positive!'
            end if
         
         end if
         ke_orb_last = ke_orbit
         ke_spin_last = ke_spin
         pe_last = pe
         Eorbit_last = Eorbit
      end associate
      return

      end subroutine io_conservation_report
end submodule s_io_conservation_report
