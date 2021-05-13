function util_solve_linear_system(A,b,n,lerr) result(x)
   !! Author: David A. Minton
   !!
   !! Solves the linear equation of the form A*x = b for x. 
   !!   A is an (n,n) arrays
   !!   x and b are (n) arrays
   !! Uses Gaussian elimination, so will have issues if system is ill-conditioned.
   !! Uses quad precision intermidiate values, so works best on small arrays.
   use swiftest
   use module_interfaces, EXCEPT_THIS_ONE => util_solve_linear_system
   implicit none
   ! Arguments
   integer(I4B),             intent(in) :: n
   real(DP), dimension(:,:), intent(in) :: A
   real(DP), dimension(:),   intent(in) :: b
   logical,                  intent(out) :: lerr
   ! Result
   real(DP), dimension(n)  :: x
   ! Internals
   real(QP), dimension(:), allocatable :: qx

   qx = solve_wbs(ge_wpp(real(A, kind=QP), real(b, kind=QP)),lerr)
   x = real(qx, kind=DP)
   return

   contains

      function solve_wbs(u, lerr) result(x) ! solve with backward substitution
         !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
         implicit none
         ! Arguments
         real(QP), intent(in), dimension(:,:), allocatable  :: u
         logical, intent(out)  :: lerr
         ! Result
         real(QP), dimension(:), allocatable :: x
         ! Internals
         integer(I4B)             :: i,n
         real(DP), parameter :: epsilon = 10 * tiny(1._DP) 

         lerr = any(abs(u(:,:)) < epsilon)
         if (lerr) return

         n = size(u,1)
         allocate(x(n))
         do i = n,1,-1 
            x(i) = (u(i, n + 1) - sum(u(i, i + 1:n) * x(i + 1:n))) / u(i, i)
         end do
         return
      end function solve_wbs

      function ge_wpp(A, b) result(u) ! gaussian eliminate with partial pivoting
         !! Solve  Ax=b  using Gaussian elimination then backwards substitution.
         !!   A being an n by n matrix.
         !!   x and b are n by 1 vectors. 
         !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
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

   end function util_solve_linear_system