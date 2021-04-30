function symba_casehitandrun (symba_plA, idx_parent, nmergeadd, mergeadd_list, name, x, v, mass, radius, L_spin, Ip,  &
                                        mass_res, param, Qloss) result(status)
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
   ! Arguments
   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   type(symba_pl), intent(inout)             :: symba_pla
   integer(I4B), dimension(:), intent(in)    :: name
   real(DP), dimension(:),   intent(in)      :: mass, radius, mass_res
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   type(user_input_parameters),intent(inout) :: param
   integer(I4B), dimension(2), intent(inout) :: idx_parent
   real(DP), intent(inout)                   :: Qloss
   ! Result
   integer(I4B)                              :: status
   ! Internals
   integer(I4B)                            :: i, nfrag, jproj, jtarg
   real(DP)                                :: mtot, avg_dens
   real(DP), dimension(NDIM)               :: xcom, vcom
   real(DP), dimension(2)                  :: vol
   real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
   real(DP), dimension(:), allocatable     :: m_frag, rad_frag
   integer(I4B), dimension(:), allocatable :: name_frag
   logical                                 :: lpure

   mtot = sum(mass(:))
   xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
   vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot
   lpure = .false.
  
   ! The largest body will stay untouched
   if (mass(1) > mass(2)) then
      jtarg = 1
      jproj = 2
   else
      jtarg = 2
      jproj = 1
   end if

   if (mass_res(2) > 0.9_DP * mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
      write(*,*) 'Pure hit and run. No new fragments generated.'
      nfrag = 0
      lpure = .true.
   else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
      nfrag = 10
      lpure = .false.
      allocate(m_frag(nfrag))
      allocate(name_frag(nfrag))
      allocate(rad_frag(nfrag))
      allocate(xb_frag(NDIM, nfrag))
      allocate(vb_frag(NDIM, nfrag))
      allocate(rot_frag(NDIM, nfrag))
      allocate(Ip_frag(NDIM, nfrag))
      m_frag(1) = mass(jtarg)
      name_frag(1) = name(jtarg)
      rad_frag(1) = radius(jtarg)
      xb_frag(:, 1) = x(:, jtarg) 
      vb_frag(:, 1) = v(:, jtarg)
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
      call symba_frag_pos(param, symba_plA, idx_parent, x, v, L_spin, Ip, mass, radius, &
                           Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, lpure, Qloss)
      if (lpure) then
         write(*,*) 'Should have been a pure hit and run instead'
         nfrag = 0
      else
         write(*,'("Generating ",I2.0," fragments")') nfrag
         do i = 1, nfrag
            name_frag(i) = param%plmaxname + i 
         end do
         param%plmaxname = name_frag(nfrag)
      end if
   end if
   if (lpure) then
      status = ACTIVE
   else
      status = HIT_AND_RUN
      ! Populate the list of new bodies
      call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  
      do i = 1, nfrag
         nmergeadd = nmergeadd + 1
         mergeadd_list%name(nmergeadd) = name_frag(i) 
         mergeadd_list%status(nmergeadd) = HIT_AND_RUN
         mergeadd_list%xb(:,nmergeadd) = xb_frag(:, i)
         mergeadd_list%vb(:,nmergeadd) = vb_frag(:, i)
         mergeadd_list%mass(nmergeadd) = m_frag(i)
         mergeadd_list%radius(nmergeadd) = rad_frag(i)
         mergeadd_list%Ip(:,nmergeadd) = Ip_frag(:, i)
         mergeadd_list%rot(:,nmergeadd) = rot_frag(:, i)
      end do 
   end if

   return 
end function symba_casehitandrun
