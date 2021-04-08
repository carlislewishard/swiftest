subroutine symba_casehitandrun (symba_plA, idx_parents, nmergeadd, mergeadd_list, name, x, v, mass, radius, L_spin, Ip, xbs, vbs, &
                                        mass_res, param)
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Create the fragments resulting from a non-catastrophic hitandrun collision
   !! 
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, EXCEPT_THIS_ONE => symba_casehitandrun
   implicit none

   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   type(symba_pl), intent(inout)             :: symba_pla
   integer(I4B), dimension(:), intent(in)    :: name
   real(DP), dimension(:),   intent(in)      :: mass, radius, xbs, vbs, mass_res
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   type(user_input_parameters),intent(inout) :: param
   integer(I4B), dimension(2), intent(inout) :: idx_parents

   integer(I4B)                            :: i, nfrag, jproj, jtarg
   real(DP)                                :: mtot, avg_dens
   real(DP), dimension(NDIM)               :: xcom, vcom
   real(DP), dimension(2)                  :: vol
   real(DP), dimension(:, :), allocatable  :: v_frag, x_frag, rot_frag, Ip_frag
   real(DP), dimension(:), allocatable     :: m_frag, rad_frag
   integer(I4B), dimension(:), allocatable :: name_frag
   logical                                 :: lmerge

   mtot = sum(mass(:))
   xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
   vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot
  
   ! The largest body will stay untouched
   if (mass(1) > mass(2)) then
      jtarg = 1
      jproj = 2
   else
      jtarg = 2
      jproj = 1
   end if

   if (mass_res(2) > 0.9_DP * mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
      nfrag = 2
      allocate(m_frag, source = mass)
      allocate(name_frag, source = name)
      allocate(rad_frag, source = radius)
      allocate(x_frag, source = x)
      allocate(v_frag, source = v)
      allocate(Ip_frag, source = Ip)
      allocate(rot_frag(NDIM, nfrag))
      do i = 1, 2
         rot_frag(:,i) = L_spin(:, i) / (Ip_frag(3, i) * m_frag(i) * rad_frag(i)**2)
      end do
      
   else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
      nfrag = 10
      allocate(m_frag(nfrag))
      allocate(name_frag(nfrag))
      allocate(rad_frag(nfrag))
      allocate(x_frag(NDIM, nfrag))
      allocate(v_frag(NDIM, nfrag))
      allocate(rot_frag(NDIM, nfrag))
      allocate(Ip_frag(NDIM, nfrag))
      m_frag(1) = mass(jtarg)
      name_frag(1) = name(jtarg)
      rad_frag(1) = radius(jtarg)
      x_frag(:, 1) = x(:, jtarg) 
      v_frag(:, 1) = v(:, jtarg)
      Ip_frag(:,1) = Ip(:, jtarg)

      ! Get mass weighted mean of Ip and average density
      vol(:) = 4._DP / 3._DP * pi * radius(:)**3
      avg_dens = mass(jproj) / vol(jproj)
      m_frag(2:nfrag) = (mtot - m_frag(1)) / (nfrag - 1) 
      rad_frag(2:nfrag) = (3 * m_frag(2:nfrag) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)
      m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))

      param%plmaxname = max(param%plmaxname, param%tpmaxname)
      
      do i = 1, nfrag
         Ip_frag(:, i) = Ip(:, jproj)
      end do

      ! Put the fragments on the circle surrounding the center of mass of the system
      call symba_frag_pos(symba_plA, idx_parents, x, v, L_spin, Ip, mass, radius, &
                           Ip_frag, m_frag, rad_frag, x_frag, v_frag, rot_frag, lmerge)
      if (lmerge) then
         write(*,*) 'Should have been a pure hit and run instead'
         nfrag = 2
         deallocate(m_frag);    allocate(m_frag, source = mass)
         deallocate(name_frag); allocate(name_frag, source = name)
         deallocate(rad_frag);  allocate(rad_frag, source = radius)
         deallocate(x_frag);    allocate(x_frag, source = x)
         deallocate(v_frag);    allocate(v_frag, source = v)
         deallocate(Ip_frag);   allocate(Ip_frag, source = Ip)
         deallocate(rot_frag);  allocate(rot_frag(NDIM, nfrag))
         do i = 1, 2
            rot_frag(:,i) = L_spin(:, i) / (Ip_frag(3, i) * m_frag(i) * rad_frag(i)**2)
         end do
      else
         do i = 1, nfrag
            name_frag(i) = param%plmaxname + i 
         end do
      end if
      param%plmaxname = name_frag(nfrag)
   end if

   ! Populate the list of new bodies
   call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  
   do i = 1, nfrag
      nmergeadd = nmergeadd + 1
      mergeadd_list%name(nmergeadd) = name_frag(i) 
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%xh(:,nmergeadd) = x_frag(:, i) - xbs(:)
      mergeadd_list%vh(:,nmergeadd) = v_frag(:, i) - vbs(:)
      mergeadd_list%mass(nmergeadd) = m_frag(i)
      mergeadd_list%radius(nmergeadd) = rad_frag(i)
      mergeadd_list%Ip(:,nmergeadd) = Ip_frag(:, i)
      mergeadd_list%rot(:,nmergeadd) = rot_frag(:, i)
   end do 

   return 
end subroutine symba_casehitandrun
