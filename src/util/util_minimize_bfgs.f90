function util_minimize_bfgs(f, N, x1, eps) result(fnum)
   !! author: David A. Minton
   !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
   !! This function implements the Broyden-Fletcher-Goldfarb-Shanno method to determine the minimum of a function of N variables.  
   !! It recieves as input:
   !!   f(x) : function of one real(DP) array variable as input
   !!   N    :  Number of variables of function f
   !!   x1   :  Initial starting value of x
   !!   eps  :  Accuracy of 1 - dimensional minimization at each step
   !! The outputs include
   !!   x    :  Final minimum (all 0 if none found)
   !! Returns
   !!   Number of function calls performed or
   !!   0 = No miniumum found
   !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
   use swiftest
   use swiftest_globals
   use module_interfaces, EXCEPT_THIS_ONE => util_minimize_bfgs
   implicit none
   ! Arguments
   integer(I4B), intent(in) :: N
   interface 
      pure function f(x) ! Objective function template
         import DP
         real(DP), dimension(:), intent(in) :: x
         real(DP) :: f
      end function 
   end interface
   real(DP), dimension(:), intent(inout) :: x1
   real(DP), intent(in) :: eps
   ! Result
   integer(I4B) :: fnum
   ! Internals
   integer(I4B) ::  i, j, k, l, conv, num
   integer(I4B), parameter :: MAXLOOP = 2000  !! Maximum number of loops before method is determined to have failed 
   real(DP), dimension(N) :: S                ! Direction vectors 
   real(DP), dimension(N) :: x0
   real(DP), dimension(N) :: Snorm           ! normalized direction 
   real(DP), dimension(N,N) :: H          ! Approximated inverse Hessian matrix 
   real(DP), dimension(N) :: grad         ! gradient of f 
   real(DP), dimension(N) :: grad0        ! old value of gradient 
   real(DP) :: astar            ! 1 - D minimized value 
   real(DP), dimension(N) :: y, P
   real(DP), dimension(N,N) :: PP, PyH, HyP
   real(DP) :: yHy, Py
   real(DP) :: xmag

   fnum = 0
   H(:,:) = reshape([((0._DP, i=1, j-1), 1._DP, (0._DP, i=j+1, N), j=1, N)], [N,N])  
   grad(:) = 0.0_DP
   do i = 1, MAXLOOP 
      xmag = norm2(x1(:))
      grad0(:) = grad(:)
      fnum = fnum + gradf(f, N, x1, grad, eps)
      if (i > 1) then
         ! set up factors for H matrix update 
         y(:) = grad(:) - grad0(:)
         P(:) = x1(:) - x0(:)
         Py = sum(P(:) * y(:))
         yHy = 0._DP
         do k = 1, N 
            yHy = yHy + y(k) * sum(H(:,k) * y(:))
         end do
         ! prevent divide by zero (convergence) 
         if (abs(Py) < tiny(Py)) return
         ! set up update 
         PyH(:,:) = 0._DP
         HyP(:,:) = 0._DP
         do k = 1, N 
            do j = 1, N
               PP(j, k) = P(j) * P(k)
               PyH(j, k) = P(j) * sum(y(:) * H(:,k))
               HyP(j, k) = P(k) * sum(y(:) * H(j,:))
            end do
         end do
         ! update H matrix 
         H(:,:) = H(:,:) + ((1._DP - yHy / Py) * PP(:, :) - PyH(:, :) - HyP(:, :)) / Py
      end if
      !check for convergence
      conv = 0
      S(:) = 0._DP
      do k = 1, N
         if (abs(grad(j)) > eps) conv = conv + 1
         S(k) = -sum(H(:,k) * grad(:))
      end do
      if (conv == 0)  return 
      ! normalize gradient 
      Snorm(:) = S(:) / norm2(S)
      num = fnum + minimize1D(f, x1, Snorm, N, eps, astar)
      if (num == fnum) then
         write(*,*) "Exiting BFGS"
         fnum = 0
         return 
      end if
      fnum = num
      ! Get new x values 
      x0(:) = x1(:)
      x1(:) = x1(:) + astar * Snorm(:)
   end do
   return 

   contains

      function gradf(f, N, x1, grad, dx) result(fnum)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! Purpose:  Estimates the gradient of a function using a central difference
         !! approximation
         !! Inputs:
         !!   f(x) :  function of a real array of x(N) as input
         !!   N    :  number of variables N
         !!   x1   :  x value array
         !!   dx   :  step size to use when calculating derivatives
         !! Output:
         !!   grad :  N sized array containing estimated gradient of f at x1
         !! Returns
         !!   Number of function calls performed
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B), intent(in) :: N
         interface 
            pure function f(x) 
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in) :: x1
         real(DP), dimension(:), intent(out) :: grad
         real(DP), intent(in) :: dx
         ! Result
         integer(I4B) :: fnum
         ! Internals
         integer(I4B) :: i, j, k
         real(DP), dimension(N) :: xp, xm

         fnum = 0
         do i = 1, N
            xp(:) = [((x1(j), j=1, k-1), x1(j) + dx, (x1(j), j=k+1, N), k=1, N)]
            xm(:) = [((x1(j), j=1, k-1), x1(j) - dx, (x1(j), j=k+1, N), k=1, N)]
            grad(i) = (f(xp) - f(xm)) / (2 * dx)
            fnum = fnum + 2
         end do
         return 
      end function gradf

      function minimize1D(f, x0, S, N, eps, astar) result(fnum)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This program find the minimum of a function of N variables in a single direction
         !! S using in sequence:
         !!    1.  A Bracketing method
         !!    2.  The golden section method
         !!    3.  A quadratic polynomial fit
         !! Inputs
         !!   f(x) : function of one real(DP) array variable as input
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   N    :  Number of variables of function f
         !!   eps  :  Accuracy of 1 - dimensional minimization at each step
         !! Output
         !!   astar      :  Final minimum along direction S
         !! Returns
         !!   Number of function calls performed or
         !!   0 = No miniumum found
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in) :: N
         interface 
            pure function f(x) ! Objective function template
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in) :: x0, S
         real(DP),               intent(in) :: eps
         real(DP),               intent(out) :: astar
         ! Result
         integer(I4B) :: fnum 
         ! Internals
         integer(I4B) :: num = 0
         real(DP), parameter :: step = 0.7_DP     !! Bracketing method step size   
         real(DP), parameter :: gam = 1.2_DP      !! Bracketing method expansion parameter   
         real(DP), parameter :: greduce = 0.2_DP  !! Golden section method reduction factor   
         real(DP), parameter :: greduce2 = 0.0_DP ! Secondary golden section method reduction factor   
         real(DP) :: alo, ahi                     !! High and low values for 1 - D minimization routines   
         real(DP), parameter :: a0 = 0.0_DP       !! Initial guess of alpha   
        
         fnum = 0
         alo = a0
         num = fnum + bracket(f, x0, S, N, alo, ahi, gam, step) 
         if (num == fnum) then
            write(*,*) "Bracketing step failed!"
            fnum = 0
            return 
         end if
         fnum = num
         if (abs(alo - ahi) < eps) then
            astar = alo
            return 
         end if
         num = fnum + golden(f, x0, S, N, alo, ahi, greduce)
         if (num == fnum) then
            write(*,*) "Golden section step failed!"
            fnum = 0
            return 
         end if
         fnum = num
         if (abs(alo - ahi) < eps) then
            astar = alo
            return 
         end if
         fnum = fnum + quadfit(f, x0, S, N, alo, ahi, eps)
         if (abs(alo - ahi) < eps) then
            astar = alo
            return 
         end if 
         ! Quadratic fit method won't converge, so finish off with another golden section   
         fnum = fnum + golden(f, x0, S, N, alo, ahi, greduce2)
         astar = (alo + ahi) / 2.0_DP
         return 
      end function minimize1D

      function n2one(f, x0, S, N, a) result(fnew)
         implicit none
         ! Arguments
         integer(I4B),           intent(in) :: N
         interface 
            pure function f(x) ! Objective function template
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in) :: x0, S
         real(DP), intent(in) :: a
         ! Return
         real(DP) :: fnew
         ! Internals
         real(DP), dimension(N) :: xnew
         integer(I4B) :: i
         
         xnew(:) = x0(:) + a * S(:)
         fnew = f(xnew(:))
         return 
      end function n2one

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      function bracket(f, x0, S, N, lo, hi, gam, step) result(fnum)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function brackets the minimum.  It recieves as input:
         !!   f(x) : function of one real(DP) array variable as input
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   lo  :  initial guess of lo bracket value
         !!   gam  :  expansion parameter
         !!   step :  step size
         !! The outputs include
         !!   lo  :  lo bracket
         !!   hi   :  hi bracket
         !! Returns
         !!   Number of function calls performed or
         !!   0 = No miniumum found
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)  :: N
         interface 
            pure function f(x) ! Objective function template
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(inout) :: lo, hi
         real(DP),               intent(in)    :: gam, step
         ! Result
         integer(I4B)                          :: fnum 
         ! Internals
         real(DP) :: a0, a1, a2, a3, da
         real(DP) :: f0, f1, f2, fcon
         integer(I4B) :: i
         integer(I4B), parameter :: MAXLOOP = 1000 ! maximum number of loops before method is determined to have failed   
         real(DP), parameter :: eps = tiny(lo) ! small number precision to test floating point equality   
         real(DP), parameter :: dela = 12.4423_DP ! arbitrary number to test if function is constant   

         ! set up initial bracket points   
         a0 =  lo
         da = step
         a1 = a0 + da
         a2 = a0 + 2 * da
         f0 = n2one(f, x0, S, N, a0)
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
         fnum = 3 ! set number off function calls   
         ! loop over bracket method until either min is bracketed method fails   
         do i = 1, MAXLOOP 
            if ((f0 > f1) .and. (f2 > f1)) then  ! minimum was found   
               lo = a0
               hi = a2
               return 
            else if ((f0 > f1) .and. (f1 < f2)) then ! function appears to decrease   
               da = da * gam
               a3 = a2 + da
               a0 = a1
               a1 = a2
               a2 = a3
               f0 = f1
               f1 = f2
               f2 = n2one(f, x0, S, N, a2)
               fnum = fnum + 1
            else if ((f0 < f1) .and. (f1 < f2)) then ! function appears to increase   
               da = da * gam
               a3 = a0 - da
               a2 = a1
               a1 = a0
               a0 = a3
               f2 = f1
               f0 = n2one(f, x0, S, N, a0)
               fnum = fnum + 1
            else if ((f0 < f1) .and. (f1 > f2)) then ! more than one minimum present, so in this case we arbitrarily choose the RHS min   
               da = da * gam
               a3 = a2 + da
               a0 = a1
               a1 = a2
               a2 = a3
               f0 = f1
               f1 = f2
               f2 = n2one(f, x0, S, N, a2)
               fnum = fnum + 1
            else if ((f0 > f1) .and. (abs(f2 - f1) <= eps)) then ! RHS equal   
               da = da * gam
               a3 = a2 + da
               a2 = a3
               f2 = n2one(f, x0, S, N, a2)
               fnum = fnum + 1
            else if ((abs(f0 - f1) < eps) .and. (f2 > f1)) then ! LHS equal   
               da = da * gam
               a3 = a0 - da
               a0 = a3
               f0 = n2one(f, x0, S, N, a0)
               fnum = fnum + 1
            else  ! all values equal stops if there is no minimum or takes RHS min if it exists   
               ! test if function itself is constant   
               fcon = n2one(f, x0, S, N, a2 + dela) !multiply be an arbitrary number to see if constant  
               fnum = fnum + 1
               if (abs(f2 - fcon) < eps) then   
                  fnum = 0
                  return ! function is constant   
               end if
               a3 = a0 + 0.5_DP * (a1 - a0)
               a0 = a1
               a1 = a2
               a2 = a3
               f0 = f1
               f1 = f2     
               a3 = a0 + 0.5_DP * (a1 - a0)
               a1 = a2
               a2 = a3
               f1 = f2
               f2 = n2one(f, x0, S, N, a2)
               fnum = fnum + 1
            end if
         end do
         fnum = 0
         return ! no minimum found   
      end function bracket

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      function golden(f, x0, S, N, lo, hi, eps) result(fnum)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses the golden section method to reduce the starting interval lo, hi by some amount sigma.  
         !! It recieves as input:
         !!   f(x) : function of one real(DP) array variable as input
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   lo   :  initial guess of lo bracket value
         !!   gam  :  expansion parameter
         !!   eps  :  reduction interval in range (0 < sigma < 1) such that:
         !!             hi(new) - lo(new) = eps * (hi(old) - lo(old))
         !! The outputs include
         !!   lo   :  lo bracket
         !!   hi   :  hi bracket
         !! Returns
         !!   Number of function calls performed or
         !!   0 = No miniumum found
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B), intent(in) :: N
         interface 
            pure function f(x) ! Objective function template
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(inout) :: lo, hi
         real(DP),               intent(in)    :: eps
         ! Result
         integer(I4B)                          :: fnum 
         ! Internals 
         real(DP), parameter :: tau = 0.5_DP * (sqrt(5.0_DP) - 1.0_DP)  ! Golden section constant   
         integer(I4B), parameter :: MAXLOOP = 40 ! maximum number of loops before method is determined to have failed (unlikely, but could occur if no minimum exists between lo and hi)   
         real(DP) :: i0 ! Initial interval value   
         real(DP) :: a1, a2
         real(DP) :: f1, f2
         integer(I4B) :: i, j

         i0 =  hi - lo
         fnum = 0
         a1 =  hi - tau * i0
         a2 =  lo + tau * i0
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
         fnum = fnum + 2
         do i = 1, MAXLOOP 
            if (abs((hi - lo) / i0) <= eps) return ! interval reduced to input amount   
            if (f2 > f1) then
               hi = a2
               a2 = a1
               f2 = f1
               a1 = hi - tau * (hi - lo)
               f1 = n2one(f, x0, S, N, a1)
               fnum = fnum + 1
            else 
               lo = a1
               a1 = a2
               f2 = f1
               a2 = hi - (1.0_DP - tau) * (hi - lo)
               f2 = n2one(f, x0, S, N, a2)
               fnum = fnum + 1
            end if
         end do
         fnum = 0
         return ! search took too many iterations - no minimum found   
      end function golden

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      function quadfit(f, x0, S, N, lo, hi, eps) result(fnum)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses a quadratic polynomial fit to  * locate the minimum of a function
         !! to some accuracy eps.  It recieves as input:
         !!   f(x) : function of one real(DP) :: variable as input
         !!   lo    :  low bracket value
         !!   hi    :  high bracket value
         !!   eps   :  desired accuracy of final minimum location
         !! The outputs include
         !!   lo   :  final minimum location
         !!   hi   :  final minimum location
         !! Returns
         !!   Number of function calls performed or
         !!   0 = No miniumum found
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B), intent(in) :: N
         interface 
            pure function f(x) ! Objective function template
               import DP
               real(DP), dimension(:), intent(in) :: x
               real(DP) :: f
            end function f
         end interface
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(inout) :: lo, hi
         real(DP),               intent(in)    :: eps
         ! Result
         integer(I4B)                          :: fnum 
         ! Internals 
         integer(I4B), parameter :: MAXLOOP = 20 ! maximum number of loops before method is determined to have failed.   
         real(DP), parameter :: small = tiny(eps) ! small number to prevent divide by zero errors   
         real(DP) :: a1, a2, a3, astar   ! three points for the polynomial fit and polynomial minimum   
         real(DP) :: f1, f2, f3, fstar   ! three function values for the polynomial and polynomial minimum   
         real(DP), dimension(3) :: row_1, row_2, row_3, rhs, soln        ! matrix for 3 equation solver (gaussian elimination)   
         real(DP), dimension(3,3) :: lhs
         real(DP) :: d1, d2, d3, aold, denom
         integer(I4B) :: i
         logical :: lerr

         ! Get initial a1, a2, a3 values   
         a1 =  lo
         a2 =  lo + 0.5_DP * (hi - lo)
         a3 =  hi
         fnum = 0
         aold = a1
         astar = a2
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
         f3 = n2one(f, x0, S, N, a3)
         fnum = fnum + 3
         do i = 1, MAXLOOP 
            ! check to see if convergence is reached and exit   
            if (abs(astar) < small) then
               denom = small 
            else 
               denom = astar
            end if
            if (abs((astar - aold) / denom) < eps) then
               lo = astar
               hi = astar
               return 
            end if
            ! Set up system for gaussian elimination equation solver   
            row_1 = [1.0_DP, a1, a1**2]
            row_2 = [1.0_DP, a2, a2**2]
            row_3 = [1.0_DP, a3, a3**2]
            rhs = [f1, f2, f3]
            lhs(1, :) = row_1
            lhs(2, :) = row_2
            lhs(3, :) = row_3
            ! Solve system of equations   
            soln(:) = util_solve_linear_system(lhs, rhs, 3, lerr)
            if (lerr) then
               write(*,*) "Could not solve polynomial on loop ", i
               write(*,'("a1 = ",f9.6," f1 = ",f9.6f)') a1, f1
               write(*,'("a2 = ",f9.6," f2 = ",f9.6f)') a2, f2
               write(*,'("a3 = ",f9.6," f3 = ",f9.6f)') a3, f3
               write(*,'("aold = ",f7.4)') aold
               fnum = 0
               return 
            end if
            aold = astar
            if (abs(soln(2)) < small) soln(2) = small
            astar =  -soln(1) / (2 * soln(2))
            fstar = n2one(f, x0, S, N, astar)
            fnum = fnum + 1
            ! keep the three closest a values to astar and discard the fourth  
            d1 = abs(a1 - astar)
            d2 = abs(a2 - astar)
            d3 = abs(a3 - astar)

            if (d1 > d2) then
               if (d1 > d3) then
                  f1 = fstar
                  a1 = astar
               else if (d3 > d2) then
                  f3 = fstar
                  a3 = astar
               end if
            else 
               if (d2 > d3) then
                  f2 = fstar
                  a2 = astar
               else if (d3 > d1) then
                  f3 = fstar
                  a3 = astar
               end if
            end if
         end do
         lo = a1
         hi = a3
         return 
      end function quadfit


end function util_minimize_bfgs