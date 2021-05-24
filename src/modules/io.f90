module io
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use swiftest_globals
   use swiftest_data_structures
   use user

   interface

      module subroutine io_getn(param,swiftest_plA,swiftest_tpA)
         type(user_input_parameters),intent(inout) :: param      !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store number of massive bodies
         type(swiftest_tp), intent(inout)  :: swiftest_tpA  !! Swiftest data structure to store number of test particles
      end subroutine io_getn

      module subroutine io_write_frame(t, swiftest_plA, swiftest_tpA, param)
         real(DP), intent(in)             :: t              !! Current time of simulation
         type(swiftest_pl), intent(inout) :: swiftest_plA   !! Swiftest massive body structure
         type(swiftest_tp), intent(inout) :: swiftest_tpA   !! Swiftest test particle structure
         type(user_input_parameters), intent(in) :: param   !! Input colleciton of user-defined parameters
      end subroutine io_write_frame

      module subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         integer(I4B), intent(in) :: iu        !! Output file unit number
         real(DP),     intent(in) :: t         !! Current time of simulation
         integer(I4B), intent(in) :: npl       !! Number of massive bodies
         integer(I4B), intent(in) :: ntp       !! Number of test particles
         integer(I4B), intent(in) :: iout_form !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in) :: out_type  !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr      
      
      module subroutine io_conservation_report(t, symba_plA, npl, j2rp2, j4rp4, param, lterminal)
         use module_symba
         real(DP),                    intent(in)    :: t            !! Current time of simulation
         type(symba_pl),              intent(inout) :: symba_plA    !! Swiftest planet data structure
         integer(I4B),                intent(in)    :: npl          !! Number of massive bodies
         real(DP),                    intent(in)    :: j2rp2, j4rp4 !! Central body oblateness terms
         type(user_input_parameters), intent(in)    :: param        !! Input colleciton of user-defined parameters
         logical,                     intent(in)    :: lterminal    !! Indicates whether to output information to the terminal screen
      end subroutine io_conservation_report

      module subroutine io_write_particle_pl(swiftest_plA, idx, param)
         implicit none
         class(swiftest_pl),          intent(in) :: swiftest_plA !! Swiftest massive body structure
         integer(I4B), dimension(:),  intent(in) :: idx       !! Array of particle indices to append to the particle file
         type(user_input_parameters), intent(in) :: param     !! Input colleciton of user-defined parameters
      end subroutine io_write_particle_pl

      module subroutine io_read_particle_pl(swiftest_plA, param)
         implicit none
         class(swiftest_pl),          intent(inout) :: swiftest_plA !! Swiftest massive body structure
         type(user_input_parameters), intent(in) :: param     !! Input colleciton of user-defined parameters
      end subroutine io_read_particle_pl
         
         
         
   end interface

end module io


