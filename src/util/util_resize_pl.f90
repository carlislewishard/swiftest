subroutine util_resize_pl(symba_plA, npl_new)
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! resize a symba massive body data structure 

! modules
   use swiftest
   use module_symba
   use module_helio
   use module_nrutil
   use module_swiftestalloc
   use module_interfaces, EXCEPT_THIS_ONE => util_resize_pl
   implicit none

! Arguments
   type(symba_pl), intent(inout) :: symba_plA
   integer(I4B), intent(in)    :: npl_new

! Internals
   type(symba_pl)          :: new_symba_plA

! Executable code

   associate(npl_old => symba_plA%helio%swiftest%nbody)
      call symba_pl_allocate(new_symba_plA, npl_new)
      if (npl_new >= npl_old) then 
         new_symba_plA%helio%swiftest%name(1:npl_old) = symba_plA%helio%swiftest%name(1:npl_old)
         new_symba_plA%helio%swiftest%status(1:npl_old) = symba_plA%helio%swiftest%status(1:npl_old)
         new_symba_plA%helio%swiftest%mass(1:npl_old) = symba_plA%helio%swiftest%mass(1:npl_old)
         new_symba_plA%helio%swiftest%radius(1:npl_old) = symba_plA%helio%swiftest%radius(1:npl_old)
         new_symba_plA%helio%swiftest%xh(:,1:npl_old) = symba_plA%helio%swiftest%xh(:,1:npl_old)
         new_symba_plA%helio%swiftest%vh(:,1:npl_old) = symba_plA%helio%swiftest%vh(:,1:npl_old)
         new_symba_plA%helio%swiftest%rhill(1:npl_old) = symba_plA%helio%swiftest%rhill(1:npl_old)
         new_symba_plA%helio%swiftest%xb(:,1:npl_old) = symba_plA%helio%swiftest%xb(:,1:npl_old)
         new_symba_plA%helio%swiftest%vb(:,1:npl_old) = symba_plA%helio%swiftest%vb(:,1:npl_old)
         new_symba_plA%helio%swiftest%rot(:,1:npl_old) = symba_plA%helio%swiftest%rot(:,1:npl_old)
         new_symba_plA%helio%swiftest%Ip(:,1:npl_old) = symba_plA%helio%swiftest%ip(:,1:npl_old)
      else
         new_symba_plA%helio%swiftest%name(1:npl_new) = symba_plA%helio%swiftest%name(1:npl_new)
         new_symba_plA%helio%swiftest%status(1:npl_new) = symba_plA%helio%swiftest%status(1:npl_new)
         new_symba_plA%helio%swiftest%mass(1:npl_new) = symba_plA%helio%swiftest%mass(1:npl_new)
         new_symba_plA%helio%swiftest%radius(1:npl_new) = symba_plA%helio%swiftest%radius(1:npl_new)
         new_symba_plA%helio%swiftest%xh(:,1:npl_new) = symba_plA%helio%swiftest%xh(:,1:npl_new)
         new_symba_plA%helio%swiftest%vh(:,1:npl_new) = symba_plA%helio%swiftest%vh(:,1:npl_new)
         new_symba_plA%helio%swiftest%rhill(1:npl_new) = symba_plA%helio%swiftest%rhill(1:npl_new)
         new_symba_plA%helio%swiftest%xb(:,1:npl_new) = symba_plA%helio%swiftest%xb(:,1:npl_new)
         new_symba_plA%helio%swiftest%vb(:,1:npl_new) = symba_plA%helio%swiftest%vb(:,1:npl_new)
         new_symba_plA%helio%swiftest%rot(:,1:npl_new) = symba_plA%helio%swiftest%rot(:,1:npl_new)
         new_symba_plA%helio%swiftest%Ip(:,1:npl_new) = symba_plA%helio%swiftest%ip(:,1:npl_new)
   
      end if
      call symba_pl_deallocate(symba_plA)
      call symba_pl_allocate(symba_plA, npl_new)
      call move_alloc(new_symba_plA%helio%swiftest%name, symba_plA%helio%swiftest%name )
      call move_alloc(new_symba_plA%helio%swiftest%status, symba_plA%helio%swiftest%status)
      call move_alloc(new_symba_plA%helio%swiftest%mass, symba_plA%helio%swiftest%mass)
      call move_alloc(new_symba_plA%helio%swiftest%radius, symba_plA%helio%swiftest%radius)
      call move_alloc(new_symba_plA%helio%swiftest%xh, symba_plA%helio%swiftest%xh)
      call move_alloc(new_symba_plA%helio%swiftest%vh, symba_plA%helio%swiftest%vh)
      call move_alloc(new_symba_plA%helio%swiftest%rhill, symba_plA%helio%swiftest%rhill)
      call move_alloc(new_symba_plA%helio%swiftest%xb , symba_plA%helio%swiftest%xb)
      call move_alloc(new_symba_plA%helio%swiftest%vb, symba_plA%helio%swiftest%vb)
      call move_alloc(new_symba_plA%helio%swiftest%rot, symba_plA%helio%swiftest%rot)
      call move_alloc(new_symba_plA%helio%swiftest%Ip, symba_plA%helio%swiftest%Ip)
      call symba_pl_deallocate(new_symba_plA)

   end associate
   return

end subroutine util_resize_pl