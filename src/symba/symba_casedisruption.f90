function symba_casedisruption (symba_plA, family, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, &
                                 mass_res, param, Qloss) result(status)
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Create the fragments resulting from a non-catastrophic disruption collision
   !! 
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, EXCEPT_THIS_ONE => symba_casedisruption
   implicit none

   ! Arguments
   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(:), intent(in)    :: family
   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   real(DP), dimension(:),   intent(inout)   :: mass, radius, mass_res
   real(DP), dimension(:,:), intent(inout)   :: x, v, L_spin, Ip
   type(user_input_parameters),intent(inout) :: param
   real(DP), intent(inout)                   :: Qloss
   ! Result
   integer(I4B)                              :: status
   ! Internals
   integer(I4B)                            :: i, istart, nfrag
   real(DP)                                :: mtot, avg_dens
   real(DP), dimension(NDIM)               :: xcom, vcom, Ip_new
   real(DP), dimension(2)                  :: vol
   real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
   real(DP), dimension(:), allocatable     :: m_frag, rad_frag
   logical                                 :: lmerge

   ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
   nfrag = 11 ! This value is set for testing. This needs to be updated such that it is calculated or set by the user
   allocate(m_frag(nfrag))
   allocate(rad_frag(nfrag))
   allocate(xb_frag(NDIM, nfrag))
   allocate(vb_frag(NDIM, nfrag))
   allocate(rot_frag(NDIM, nfrag))
   allocate(Ip_frag(NDIM, nfrag))

   mtot = sum(mass(:))
   xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
   vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

   ! Get mass weighted mean of Ip and average density
   Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mtot
   vol(:) = 4._DP / 3._DP * PI * radius(:)**3
   avg_dens = mtot / sum(vol(:))

   ! Distribute the mass among fragments, with a branch to check for the size of the second largest fragment
   m_frag(1) = mass_res(1)
   if (mass_res(2) > mass_res(1) / 3._DP) then
      m_frag(2) = mass_res(2)
      istart = 3
   else
      istart = 2
   end if
   ! Distribute remaining mass among the remaining bodies
   do i = istart, nfrag
      m_frag(i) = (mtot - sum(m_frag(1:istart - 1))) / (nfrag - istart + 1) 
   end do

   ! Distribute any residual mass if there is any and set the radius
   m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
   rad_frag(:) = (3 * m_frag(:) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)

   do i = 1, nfrag
      Ip_frag(:, i) = Ip_new(:)
   end do

   ! Put the fragments on the circle surrounding the center of mass of the system
   call symba_frag_pos(param, symba_plA, family, x, v, L_spin, Ip, mass, radius, &
                       nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lmerge)

   if (lmerge) then
      write(*,*) 'Should have been a merge instead.'
      status = symba_casemerge (symba_plA, family, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, param)
   else
      ! Populate the list of new bodies
      write(*,'("Generating ",I2.0," fragments")') nfrag
      call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  
      status = DISRUPTION
      do i = 1, nfrag
         nmergeadd = nmergeadd + 1
         symba_plA%helio%swiftest%maxid = symba_plA%helio%swiftest%maxid + 1
         mergeadd_list%id(nmergeadd) = symba_plA%helio%swiftest%maxid
         mergeadd_list%status(nmergeadd) = DISRUPTION
         mergeadd_list%xb(:,nmergeadd) = xb_frag(:, i) 
         mergeadd_list%vb(:,nmergeadd) = vb_frag(:, i) 
         mergeadd_list%mass(nmergeadd) = m_frag(i)
         mergeadd_list%radius(nmergeadd) = rad_frag(i)
         mergeadd_list%Ip(:,nmergeadd) = Ip_frag(:, i)
         mergeadd_list%rot(:,nmergeadd) = rot_frag(:, i)
         mergeadd_list%info%origin_type = "Disruption"
      end do 
   end if

   return 
end function symba_casedisruption
