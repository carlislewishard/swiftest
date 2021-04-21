submodule (io) s_io_conservation_report
contains
   module procedure io_conservation_report
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Reports the current state of energy, mass, and angular momentum conservation in a run
   use module_interfaces
   implicit none

      real(DP), dimension(NDIM), save :: Ltot_orig
      real(DP), save                  :: Eorbit_orig, Mtot_orig, Lmag_orig
      real(DP)                        :: ke, pe, Eorbit
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


      if (lfirst) then
         if (param%out_stat == "OLD") then
            open(unit = egyiu, file = energy_file, form = "formatted", status = "old", action = "write", position = "append")
         else 
            open(unit = egyiu, file = energy_file, form = "formatted", status = "replace", action = "write")
            write(egyiu,egyheader)
         end if
      end if
      call symba_energy_eucl(npl, symba_plA, j2rp2, j4rp4, ke, pe, Eorbit, Ltot_now, Mtot_now)
      Mtot_now = symba_plA%helio%swiftest%Mescape + Mtot_now
      Ltot_now(:) = symba_plA%helio%swiftest%Lescape(:) + Ltot_now(:)
      if (lfirst) then
         Eorbit_orig = Eorbit
         Mtot_orig = Mtot_now
         Ltot_orig(:) = Ltot_now(:)
         Lmag_orig = norm2(Ltot_orig(:))
         lfirst = .false.
      end if

      write(egyiu,egyfmt) t, Eorbit, symba_plA%helio%swiftest%Ecollisions, Ltot_now, Mtot_now
      flush(egyiu)
      if (.not.lfirst .and. lterminal) then 
         Lmag_now = norm2(Ltot_now)
         Lerror = (Lmag_now - Lmag_orig) / Lmag_orig
         Eorbit_error = (Eorbit - Eorbit_orig) / abs(Eorbit_orig)
         Ecoll_error = -symba_plA%helio%swiftest%Ecollisions / abs(Eorbit_orig)
         Etotal_error = (Eorbit - (Eorbit_orig - symba_plA%helio%swiftest%Ecollisions)) / abs(Eorbit_orig)
         Merror = (Mtot_now - Mtot_orig) / Mtot_orig
         write(*, egytermfmt) Lerror, Ecoll_error, Etotal_error, Merror
      end if

   end procedure io_conservation_report
end submodule s_io_conservation_report
