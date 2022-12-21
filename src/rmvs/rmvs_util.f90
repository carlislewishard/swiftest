!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(rmvs) s_rmvs_util
   use swiftest
contains

   module subroutine rmvs_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(rmvs_pl),                  intent(inout) :: self         !! RMVS massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (rmvs_pl)
         associate(nold => self%nbody, nsrc => source%nbody)
            call swiftest_util_append(self%nenc, source%nenc, nold, nsrc, lsource_mask)
            call swiftest_util_append(self%tpenc1P, source%tpenc1P, nold, nsrc, lsource_mask)
            call swiftest_util_append(self%plind, source%plind, nold, nsrc, lsource_mask)

            ! The following are not implemented as RMVS doesn't make use of fill operations on pl type
            ! So they are here as a placeholder in case someone wants to extend the RMVS class for some reason
            !call swiftest_util_append(self%outer, source%outer, nold, nsrc, lsource_mask)
            !call swiftest_util_append(self%inner, source%inner, nold, nsrc, lsource_mask)
            !call swiftest_util_append(self%planetocentric, source%planetocentric, nold, nsrc, lsource_mask)

            call whm_util_append_pl(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class rmvs_pl or its descendents!"
         call swiftest_util_exit(FAILURE)
      end select

      return
   end subroutine rmvs_util_append_pl


   module subroutine rmvs_util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from test particle object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(rmvs_tp),                  intent(inout) :: self         !! RMVS test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (rmvs_tp)
         associate(nold => self%nbody, nsrc => source%nbody)
            call swiftest_util_append(self%lperi, source%lperi, nold, nsrc, lsource_mask)
            call swiftest_util_append(self%plperP, source%plperP, nold, nsrc, lsource_mask)
            call swiftest_util_append(self%plencP, source%plencP, nold, nsrc, lsource_mask)

            call swiftest_util_append_tp(self, source, lsource_mask)  ! Note: whm_tp does not have its own append method, so we skip back to the base class
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class rmvs_tp or its descendents!"
         call swiftest_util_exit(FAILURE)
      end select

      return
   end subroutine rmvs_util_append_tp


   module subroutine rmvs_util_dealloc_cb(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Argument
      class(rmvs_cb),  intent(inout) :: self !! RMVS central body object

      if (allocated(self%outer)) deallocate(self%outer)
      if (allocated(self%inner)) deallocate(self%inner)

      return
   end subroutine rmvs_util_dealloc_cb


   module subroutine rmvs_util_dealloc_interp(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Argument
      class(rmvs_interp),  intent(inout) :: self !! RMVS interpolated system variables object
      
      if (allocated(self%x)) deallocate(self%x)
      if (allocated(self%v)) deallocate(self%v)
      if (allocated(self%aobl)) deallocate(self%aobl)
      if (allocated(self%atide)) deallocate(self%atide)

      return
   end subroutine rmvs_util_dealloc_interp


   module subroutine rmvs_util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Argumente
      class(rmvs_pl),  intent(inout) :: self !! RMVS massive body object

      if (allocated(self%outer)) deallocate(self%outer)
      if (allocated(self%inner)) deallocate(self%inner)
      if (allocated(self%nenc))  deallocate(self%nenc)
      if (allocated(self%planetocentric)) deallocate(self%planetocentric)

      call whm_util_dealloc_pl(self)

      return
   end subroutine rmvs_util_dealloc_pl


   module subroutine rmvs_util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Argument
      class(rmvs_tp),  intent(inout) :: self !! RMVS test particle object

      if (allocated(self%lperi)) deallocate(self%lperi)
      if (allocated(self%plperP)) deallocate(self%plperP)
      if (allocated(self%plencP)) deallocate(self%plencP)
      if (allocated(self%rheliocentric)) deallocate(self%rheliocentric)
      call self%cb_heliocentric%dealloc()

      call swiftest_util_dealloc_tp(self)

      return
   end subroutine rmvs_util_dealloc_tp


   module subroutine rmvs_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new RMVS massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (rmvs_pl)
            call swiftest_util_fill(keeps%nenc, inserts%nenc, lfill_list)
            call swiftest_util_fill(keeps%tpenc1P, inserts%tpenc1P, lfill_list)
            call swiftest_util_fill(keeps%plind, inserts%plind, lfill_list)

            ! The following are not implemented as RMVS doesn't make use of fill operations on pl type
            ! So they are here as a placeholder in case someone wants to extend the RMVS class for some reason
            !call swiftest_util_fill(keeps%outer, inserts%outer, lfill_list)
            !call swiftest_util_fill(keeps%inner, inserts%inner, lfill_list)
            !call swiftest_util_fill(keeps%planetocentric, inserts%planetocentric, lfill_list)

            call whm_util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class rmvs_pl or its descendents!"
            call swiftest_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine rmvs_util_fill_pl


   module subroutine rmvs_final_cb(self)
      !! author: David A. Minton
      !!
      !! Finalize the RMVS massive body object - deallocates all allocatables
      implicit none
      ! Arguments
      type(rmvs_cb),  intent(inout) :: self !! RMVS central body object

      call self%dealloc()

      return
   end subroutine rmvs_final_cb


   module subroutine rmvs_final_interp(self)
      !! author: David A. Minton
      !!
      !! Finalize the RMVS nbody system object - deallocates all allocatables
      implicit none
      ! Arguments
      type(rmvs_interp),  intent(inout) :: self !! RMVS nbody system object

      call self%dealloc()

      return
   end subroutine rmvs_final_interp


   module subroutine rmvs_final_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the RMVS massive body object - deallocates all allocatables
      implicit none
      ! Arguments
      type(rmvs_pl),  intent(inout) :: self !! RMVS massive body object

      call self%dealloc()

      return
   end subroutine rmvs_final_pl


   module subroutine rmvs_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the RMVS nbody system object - deallocates all allocatables
      implicit none
      ! Arguments
      type(rmvs_nbody_system),  intent(inout) :: self !! RMVS nbody system object

      if (allocated(self%vbeg)) deallocate(self%vbeg)
      call whm_final_system(self%whm_nbody_system)

      return
   end subroutine rmvs_final_system


   module subroutine rmvs_final_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the RMVS test particle object - deallocates all allocatables
      implicit none
      ! Arguments
      type(rmvs_tp),  intent(inout) :: self !! RMVS test particle object

      call self%dealloc()

      return
   end subroutine rmvs_final_tp


   module subroutine rmvs_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new RMVS test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(rmvs_tp),        intent(inout) :: self       !! RMVS test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (rmvs_tp)
            call swiftest_util_fill(keeps%lperi, inserts%lperi, lfill_list)
            call swiftest_util_fill(keeps%plperP, inserts%plperP, lfill_list)
            call swiftest_util_fill(keeps%plencP, inserts%plencP, lfill_list)
            
            call swiftest_util_fill_tp(keeps, inserts, lfill_list) ! Note: whm_tp does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class rmvs_tp or its descendents!"
            call swiftest_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine rmvs_util_fill_tp


   module subroutine rmvs_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a massive body object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout) :: self  !! RMVS massive body object
      integer(I4B),   intent(in)    :: nnew  !! New size neded

      call swiftest_util_resize(self%nenc, nnew)
      call swiftest_util_resize(self%tpenc1P, nnew)
      call swiftest_util_resize(self%plind, nnew)

      ! The following are not implemented as RMVS doesn't make use of resize operations on pl type
      ! So they are here as a placeholder in case someone wants to extend the RMVS class for some reason
      !call swiftest_util_resize(self%outer, nnew)
      !call swiftest_util_resize(self%inner, nnew)
      !call swiftest_util_resize(self%planetocentric, nnew)

      call whm_util_resize_pl(self, nnew)
      return
   end subroutine rmvs_util_resize_pl


   module subroutine rmvs_util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a test particle object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(rmvs_tp), intent(inout) :: self  !! RMVS test particle object
      integer(I4B),   intent(in)    :: nnew  !! New size neded

      call swiftest_util_resize(self%lperi, nnew)
      call swiftest_util_resize(self%plperP, nnew)
      call swiftest_util_resize(self%plencP, nnew)
      call swiftest_util_resize(self%rheliocentric, nnew)

      call swiftest_util_resize_tp(self, nnew)

      return
   end subroutine rmvs_util_resize_tp


   module subroutine rmvs_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a RMVS massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout) :: self       !! RMVS massive body object
      character(*),   intent(in)     :: sortby    !! Sorting attribute
      logical,        intent(in)     :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B) :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("nenc")
            call swiftest_util_sort(direction * pl%nenc(1:npl), ind)
         case("tpenc1P")
            call swiftest_util_sort(direction * pl%tpenc1P(1:npl), ind)
         case("plind")
            call swiftest_util_sort(direction * pl%plind(1:npl), ind)
         case("outer", "inner", "planetocentric", "lplanetocentric")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call whm_util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate
      return
   end subroutine rmvs_util_sort_pl


   module subroutine rmvs_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a RMVS test particle object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(rmvs_tp), intent(inout) :: self      !! RMVS test particle object
      character(*),   intent(in)    :: sortby    !! Sorting attribute
      logical,        intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                            :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("plperP")
            call swiftest_util_sort(direction * tp%plperP(1:ntp), ind)
         case("plencP")
            call swiftest_util_sort(direction * tp%plencP(1:ntp), ind)
         case("lperi", "cb_heliocentric", "rheliocentric", "index", "ipleP", "lplanetocentric")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class (*NOTE whm_tp does not need its own sort method, so we go straight to the swiftest_tp method)
            call swiftest_util_sort_tp(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)

      end associate
      return
   end subroutine rmvs_util_sort_tp

   module subroutine rmvs_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange RMVS massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(rmvs_pl),               intent(inout) :: self !! RMVS massive body object
      integer(I4B),   dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         call swiftest_util_sort_rearrange(pl%nenc, ind, npl)
         call swiftest_util_sort_rearrange(pl%tpenc1P, ind, npl)
         call swiftest_util_sort_rearrange(pl%plind, ind, npl)
         call swiftest_util_sort_rearrange_pl(pl,ind)
      end associate

      return
   end subroutine rmvs_util_sort_rearrange_pl


   module subroutine rmvs_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange RMVS test particle object in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self !! RMVS test particle object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         call swiftest_util_sort_rearrange(tp%lperi, ind, ntp)
         call swiftest_util_sort_rearrange(tp%plperP, ind, ntp)
         call swiftest_util_sort_rearrange(tp%plencP, ind, ntp)
         call swiftest_util_sort_rearrange(tp%rheliocentric, ind, ntp)
         call swiftest_util_sort_rearrange_tp(tp,ind)
      end associate

      return
   end subroutine rmvs_util_sort_rearrange_tp
   

   module subroutine rmvs_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine discard_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_pl),        intent(inout) :: self         !! RMVS massive body body object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not

      associate(keeps => self)
         select type(discards)
         class is (rmvs_pl)
            call swiftest_util_spill(keeps%nenc, discards%nenc, lspill_list, ldestructive)
            call swiftest_util_spill(keeps%tpenc1P, discards%tpenc1P, lspill_list, ldestructive)
            call swiftest_util_spill(keeps%plind, discards%plind, lspill_list, ldestructive)

            call whm_util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class rmvs_pl or its descendents!"
            call swiftest_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine rmvs_util_spill_pl

   
   module subroutine rmvs_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_tp),        intent(inout) :: self        !! RMVS test particle object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not

      associate(keeps => self)
         select type(discards)
         class is (rmvs_tp)
            call swiftest_util_spill(keeps%lperi, discards%lperi, lspill_list, ldestructive)
            call swiftest_util_spill(keeps%plperP, discards%plperP, lspill_list, ldestructive)
            call swiftest_util_spill(keeps%plencP, discards%plencP, lspill_list, ldestructive)

            call swiftest_util_spill_tp(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class rmvs_tp or its descendents!"
            call swiftest_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine rmvs_util_spill_tp

end submodule s_rmvs_util
