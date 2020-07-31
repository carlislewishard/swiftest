function symba_mergeadd_eoffset(npl, symba_plA, mergeadd_list, mergesub_list, addi, addf, subi, subf, param) result(eoffset)
   !! author: David A. Minton
   !!
   !! Calculates the difference in total energy of the system from bodies added/removed
   !! during a collision
   !!  
   use swiftest
   USE module_symba
   use module_interfaces, except_this_one => symba_mergeadd_eoffset
   implicit none

   integer(I4B), intent(in)                :: npl
   type(symba_pl), intent(in)              :: symba_plA
   type(symba_merger), intent(in)          :: mergeadd_list, mergesub_list
   integer(I4B), intent(in)                :: addi, addf, subi, subf
   type(user_input_parameters), intent(in) :: param
   real(DP)                                :: eoffset

   integer(I4B)              :: i, j
   real(DP)                  :: ke, pe, mass, rmag, v2, oblpot, rinv2, t0, t1, t2, t3, p2, p4, mu, irh
   real(DP), dimension(NDIM) :: dx, x, v, vbs
   logical  :: keeper


   ke = 0.0_DP
   vbs(:) = symba_plA%helio%swiftest%vb(:, 1)

   ! Remove the KE of the subtracted bodies
   !$omp simd private(x,v,v2,mass) reduction(-:ke)
   do i = subi, subf 
      x(:) = mergesub_list%xh(:, i)
      v(:) = mergesub_list%vh(:, i) + vbs(:)
      mass = mergesub_list%mass(i)
      v2 = dot_product(v(:), v(:))
      ke = ke - 0.5_DP * mass * v2
   end do

   !Add the KE of the added bodies
   !!$omp simd private(x,v,v2,mass) reduction(+:ke)
   do i = addi, addf 
      x(:) = mergeadd_list%xh(:, i)
      v(:) = mergeadd_list%vh(:, i) + vbs(:)
      mass = mergeadd_list%mass(i)
      v2 = dot_product(v(:), v(:))
      ke = ke + 0.5_DP * mass * v2
   end do

   pe = 0.0_DP
   !!$omp parallel do default(private) &
   !!$omp shared (symba_plA, npl, mergesub_list, mergeadd_list, subi, subf, addi, addf) &
   !!$omp reduction (+:pe)
   do i = 1, npl
      keeper = .true.
      do j = subi, subf
         if (mergesub_list%index_ps(j) /= i) then ! Don't do any calculations if the body is 
                                                  ! the same as the one being subtracted
            dx(:) = mergesub_list%xh(:, j) - symba_plA%helio%swiftest%xh(:, i) 
            rmag = norm2(dx(:)) 
            if (rmag > tiny(rmag)) pe = pe + symba_plA%helio%swiftest%mass(i) * mergesub_list%mass(j) / rmag
         else
            keeper = .false. ! This is a body being removed
         end if
      end do
      if (keeper) then
         do j = addi, addf
            dx(:) = mergeadd_list%xh(:, j) - symba_plA%helio%swiftest%xh(:, i) 
            rmag = norm2(dx(:)) 
            if (rmag > tiny(rmag)) pe = pe - symba_plA%helio%swiftest%mass(i) * mergeadd_list%mass(j) / rmag
         end do
      end if
   end do
   !!$omp end parallel do

   if (param%j2rp2 /= 0.0_DP) then
      !!$omp simd private(rmag)
      oblpot = 0.0_DP
      mu = symba_plA%helio%swiftest%mass(1)
      do i = subi, subf
         rmag = norm2(mergesub_list%xh(:,i))
         irh = 1.0_DP / rmag
         rinv2 = irh**2
         t0 = mu * mergesub_list%mass(i) * rinv2 * irh
         t1 = param%j2rp2
         t2 = mergesub_list%xh(3, i)**2 * rinv2
         t3 = param%j4rp4 * rinv2
         p2 = 0.5_DP * (3 * t2 - 1.0_DP)
         p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
         oblpot = oblpot - t0 * (t1 * p2 + t3 * p4)
      end do

      do i = addi, addf
         rmag = norm2(mergeadd_list%xh(:,i))
         irh = 1.0_DP / rmag
         rinv2 = irh**2
         t0 = mu * mergeadd_list%mass(i) * rinv2 * irh
         t1 = param%j2rp2
         t2 = mergeadd_list%xh(3, i)**2 * rinv2
         t3 = param%j4rp4 * rinv2
         p2 = 0.5_DP * (3 * t2 - 1.0_DP)
         p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
         oblpot = oblpot + t0 * (t1 * p2 + t3 * p4)
      end do 
      pe = pe + oblpot
   end if

   eoffset = ke + pe

   return

end function symba_mergeadd_eoffset