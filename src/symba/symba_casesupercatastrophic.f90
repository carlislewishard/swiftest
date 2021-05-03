function symba_casesupercatastrophic (symba_plA, family, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, &
                                        mass_res, param, Qloss) result(status)
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Create the fragments resulting from a supercatastrophic collision
   !! 
   use swiftest
   use module_helio
   use module_symba
   use module_swiftestalloc 
   use module_interfaces, EXCEPT_THIS_ONE => symba_casesupercatastrophic
   implicit none
   ! Arguments
   type(symba_pl), intent(inout)             :: symba_plA
   integer(I4B), dimension(:), intent(in)    :: family
   integer(I4B), intent(inout)               :: nmergeadd
   type(symba_merger), intent(inout)         :: mergeadd_list
   real(DP), dimension(:,:), intent(in)      :: x, v, L_spin, Ip
   real(DP), dimension(:),   intent(in)      :: mass, radius, mass_res
   type(user_input_parameters),intent(inout) :: param
   real(DP), intent(inout)                   :: Qloss
   ! Result
   integer(I4B)                              :: status
   ! Internals
   integer(I4B)                            :: i, idstart, nfrag
   real(DP)                                :: mtot, avg_dens, min_frag_mass
   real(DP), dimension(NDIM)               :: xcom, vcom
   real(DP), dimension(2)                  :: vol
   real(DP), dimension(NDIM)               :: Ip_new
   real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
   real(DP), dimension(:), allocatable     :: m_frag, rad_frag
   logical                                 :: lmerge
  
   ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
   nfrag = 20 ! This value is set for testing. This needs to be updated such that it is calculated or set by the user
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
   vol(:) = 4._DP / 3._DP * pi * radius(:)**3
   avg_dens = mtot / sum(vol(:))

   ! If we are adding the first and largest fragment (lr), check to see if its mass is SMALLER than an equal distribution of 
   ! mass between all fragments. If so, we will just distribute the mass equally between the fragments
   min_frag_mass = mtot / nfrag
   if (mass_res(1) < min_frag_mass) then
      m_frag(:) = min_frag_mass
   else
      m_frag(1) = mass_res(1)
      m_frag(2:nfrag) = (mtot - mass_res(1)) / (nfrag - 1)
   end if
   ! Distribute any residual mass if there is any and set the radius
   m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
   rad_frag(:) = (3 * m_frag(:) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)

   do i = 1, nfrag
      Ip_frag(:, i) = Ip_new(:)
   end do

   ! Put the fragments on the circle surrounding the center of mass of the system
   call symba_frag_pos(param, symba_plA, family, x, v, L_spin, Ip, mass, radius, &
                        Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, lmerge, Qloss)

   if (lmerge) then
      write(*,*) 'Should have been a merge instead.'
      status = symba_casemerge (symba_plA, family, nmergeadd, mergeadd_list, x, v, mass, radius, L_spin, Ip, param)
   else
      write(*,'("Generating ",I2.0," fragments")') nfrag
      status = SUPERCATASTROPHIC
      ! Populate the list of new bodies
      call symba_merger_size_check(mergeadd_list, nmergeadd + nfrag)  
      associate(npl => symba_plA%helio%swiftest%nbody)
         idstart = maxval([symba_plA%helio%swiftest%id(1:npl), mergeadd_list%id(1:nmergeadd)])
      end associate
      do i = 1, nfrag
         nmergeadd = nmergeadd + 1
         mergeadd_list%id(nmergeadd) = idstart + i
         mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
         mergeadd_list%xb(:,nmergeadd) = xb_frag(:, i) 
         mergeadd_list%vb(:,nmergeadd) = vb_frag(:, i)
         mergeadd_list%mass(nmergeadd) = m_frag(i)
         mergeadd_list%radius(nmergeadd) = rad_frag(i)
         mergeadd_list%Ip(:,nmergeadd) = Ip_frag(:, i)
         mergeadd_list%rot(:,nmergeadd) = rot_frag(:, i)
      end do 
   end if
   
   return 
end function symba_casesupercatastrophic
