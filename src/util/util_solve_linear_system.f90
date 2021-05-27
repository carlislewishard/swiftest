function util_solve_linear_system_d(A,b,n,lerr) result(x)
   !! Author: David A. Minton
   !!
   !! Solves the linear equation of the form A*x = b for x. 
   !!   A is an (n,n) arrays
   !!   x and b are (n) arrays
   !! Uses Gaussian elimination, so will have issues if system is ill-conditioned.
   !! Uses quad precision intermidiate values, so works best on small arrays.
   use, intrinsic :: ieee_exceptions
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => util_solve_linear_system_d
   implicit none
   ! Arguments
   integer(I4B),             intent(in)  :: n
   real(DP), dimension(:,:), intent(in)  :: A
   real(DP), dimension(:),   intent(in)  :: b
   logical,                  intent(out) :: lerr
   ! Result
   real(DP), dimension(n)                :: x
   ! Internals
   real(QP), dimension(:), allocatable :: qx
   type(ieee_status_type) :: original_fpe_status
   logical, dimension(:), allocatable :: fpe_flag 

   call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
   call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
   allocate(fpe_flag(size(ieee_usual)))

   qx = solve_wbs(ge_wpp(real(A, kind=QP), real(b, kind=QP)))

   call ieee_get_flag(ieee_usual, fpe_flag)
   lerr = any(fpe_flag) 
   if (lerr) then
      x = 0.0_DP
      write(*,*) 'fpe in util_solve_linear_system'
   else
      x = real(qx, kind=DP)
   end if
   call ieee_set_status(original_fpe_status)

   return
end function util_solve_linear_system_d

function util_solve_linear_system_q(A,b,n,lerr) result(x)
   !! Author: David A. Minton
   !!
   !! Solves the linear equation of the form A*x = b for x. 
   !!   A is an (n,n) arrays
   !!   x and b are (n) arrays
   !! Uses Gaussian elimination, so will have issues if system is ill-conditioned.
   !! Uses quad precision intermidiate values, so works best on small arrays.
   use, intrinsic :: ieee_exceptions
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => util_solve_linear_system_q
   implicit none
   ! Arguments
   integer(I4B),             intent(in) :: n
   real(QP), dimension(:,:), intent(in) :: A
   real(QP), dimension(:),   intent(in) :: b
   logical,                  intent(out) :: lerr
   ! Result
   real(QP), dimension(n)  :: x
   type(ieee_status_type) :: original_fpe_status
   logical, dimension(:), allocatable :: fpe_flag 

   call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
   call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
   allocate(fpe_flag(size(ieee_usual)))

   x = solve_wbs(ge_wpp(A, b))

   call ieee_get_flag(ieee_usual, fpe_flag)
   lerr = any(fpe_flag) 
   if (lerr) x = 0.0_DP
   call ieee_set_status(original_fpe_status) 

   return
end function util_solve_linear_system_q

function solve_wbs(u) result(x) ! solve with backward substitution
   !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
   use, intrinsic :: ieee_exceptions
   use swiftest
   implicit none
   ! Arguments
   real(QP), intent(in), dimension(:,:), allocatable  :: u
   ! Result
   real(QP), dimension(:), allocatable :: x
   ! Internals
   integer(I4B)             :: i,n

   n = size(u, 1)
   if (allocated(x)) deallocate(x)
   if (.not.allocated(x)) allocate(x(n))
   if (any(abs(u) < tiny(1._DP))) then 
      x(:) = 0._DP
      return
   end if
   forall (i=n:1:-1) x(i) = (u(i, n + 1) - sum(u(i, i + 1:n) * x(i + 1:n))) / u(i, i)
   return
end function solve_wbs

function ge_wpp(A, b) result(u) ! gaussian eliminate with partial pivoting
   !! Solve  Ax=b  using Gaussian elimination then backwards substitution.
   !!   A being an n by n matrix.
   !!   x and b are n by 1 vectors. 
   !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
   use, intrinsic :: ieee_exceptions
   use swiftest
   implicit none
   ! Arguments
   real(QP), dimension(:,:), intent(in) :: A
   real(QP), dimension(:),   intent(in) :: b
   ! Result
   real(QP), dimension(:,:), allocatable :: u
   ! Internals
   integer(I4B) :: i,j,n,p
   real(QP)     ::  upi

   n = size(a, 1)
   allocate(u(n, (n + 1)))
   u = reshape([A, b], [n, n + 1])
   call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
   do j = 1, n
      p = maxloc(abs(u(j:n, j)), 1) + j - 1 ! maxloc returns indices between (1, n - j + 1)
      if (p /= j) u([p, j], j) = u([j, p], j)
      u(j + 1:, j) = u(j + 1:, j) / u(j, j)
      do i = j + 1, n + 1
         upi = u(p, i)
         if (p /= j) u([p, j], i) = u([j, p], i)
         u(j + 1:n, i) = u(j + 1:n, i) - upi * u(j + 1:n, j)
      end do
   end do
   return
end function ge_wpp
