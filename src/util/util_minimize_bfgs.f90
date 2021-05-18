function util_minimize_bfgs(f, N, x0_d, eps_d, lerr) result(x1_d)
   !! author: David A. Minton
   !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
   !! This function implements the Broyden-Fletcher-Goldfarb-Shanno method to determine the minimum of a function of N variables.  
   !! It recieves as input:
   !!   f%eval(x) : lambda function object containing the objective function as the eval metho
   !!   N    :  Number of variables of function f
   !!   x0   :  Initial starting value of x
   !!   eps  :  Accuracy of 1 - dimensional minimization at each step
   !! The outputs include
   !!   lerr :  Returns .true. if it could not find the minimum
   !! Returns
   !!   x1   :  Final minimum (all 0 if none found)
   !!   0 = No miniumum found
   !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
   use, intrinsic :: ieee_exceptions
   use swiftest
   use swiftest_globals
   use module_interfaces, EXCEPT_THIS_ONE => util_minimize_bfgs
   implicit none
   ! Arguments
   integer(I4B),           intent(in)    :: N
   class(lambda_obj),      intent(in)    :: f
   real(DP), dimension(:), intent(in)    :: x0_d
   real(DP),               intent(in)    :: eps_d
   logical,                intent(out)   :: lerr
   ! Result
   real(DP), dimension(:), allocatable :: x1_d
   ! Internals
   integer(I4B) ::  i, j, k, l, conv, num
   integer(I4B), parameter :: MAXLOOP = 500 !! Maximum number of loops before method is determined to have failed 
   real(QP), parameter     :: gradeps = 1e-8_QP !! Tolerance for gradient calculations
   real(QP), dimension(N) :: S               !! Direction vectors 
   real(QP), dimension(N) :: Snorm           !! normalized direction 
   real(QP), dimension(N,N) :: H             !! Approximated inverse Hessian matrix 
   real(QP), dimension(N) :: grad1           !! gradient of f 
   real(QP), dimension(N) :: grad0           !! old value of gradient 
   real(QP) :: astar                         !! 1D minimized value 
   real(QP), dimension(N) :: y, P, x0, x1
   real(QP), dimension(N,N) :: PP, PyH, HyP
   real(QP) :: yHy, Py, eps
   type(ieee_status_type) :: original_fpe_status
   logical, dimension(:), allocatable :: fpe_flag 

   call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
   call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
   allocate(fpe_flag(size(ieee_usual)))

   lerr = .false.
   eps = real(eps_d, kind=QP)
   x0 = real(x0_d, kind=QP)
   allocate(x1_d, source=x0_d)
   x1 = x0_d
   ! Initialize approximate Hessian with the identity matrix (i.e. begin with method of steepest descent) 
   ! Get initial gradient and initialize arrays for updated values of gradient and x
   H(:,:) = reshape([((0._QP, i=1, j-1), 1._QP, (0._QP, i=j+1, N), j=1, N)], [N,N])  
   grad0 = gradf(f, N, x0(:), gradeps)
   grad1(:) = grad0(:)
   do i = 1, MAXLOOP 
      !check for convergence
      conv = 0
      S(:) = 0._QP
      do k = 1, N
         if (abs(grad1(k)) > eps) conv = conv + 1
         S(k) = -sum(H(:,k) * grad1(:))
      end do
      if (conv == 0) then
         write(*,*) "BFGS converged on gradient after ",i," iterations" 
         exit 
      end if
      ! normalize gradient 
      Snorm(:) = S(:) / norm2(S)
      astar = minimize1D(f, x1, Snorm, N, gradeps, lerr)
      if (lerr) then
         write(*,*) "Exiting BFGS with error in minimize1D step"
         exit
      end if
      ! Get new x values 
      P(:) = astar * Snorm(:) 
      x1(:) = x1(:) + P(:)
      ! Calculate new gradient
      grad0(:) = grad1(:)
      grad1 = gradf(f, N, x1, gradeps)
      y(:) = grad1(:) - grad0(:)
      Py = sum(P(:) * y(:))
      ! set up factors for H matrix update 
      yHy = 0._QP
      do k = 1, N 
         do j = 1, N
            yHy = yHy + y(j) * H(j,k) * y(k)
         end do
      end do
      ! prevent divide by zero (convergence) 
      if (abs(Py) < tiny(Py)) then
         write(*,*) "BFGS Converged on tiny Py after ",i," iterations"
         exit
      end if
      ! set up update 
      PyH(:,:) = 0._QP
      HyP(:,:) = 0._QP
      do k = 1, N 
         do j = 1, N
            PP(j, k) = P(j) * P(k)
            do l = 1, N
               PyH(j, k) = PyH(j, k) + P(j) * y(l) * H(l,k)
               HyP(j, k) = HyP(j, k) + P(k) * y(l) * H(j,l)
            end do
         end do
      end do
      ! update H matrix 
      H(:,:) = H(:,:) + ((1._QP - yHy / Py) * PP(:,:) - PyH(:,:) - HyP(:,:)) / Py
      if (any(H(:,:) > sqrt(huge(1._QP)) / N)) then
         write(*,*) 'BFGS did not converge after ',i,'iterations: H too big'
         exit
      end if
      ! Stop everything if there are any exceptions to allow the routine to fail gracefully
      call ieee_get_flag(ieee_usual, fpe_flag)
      if (any(fpe_flag)) exit 
      if (i == MAXLOOP) write(*,*) "BFGS ran out of loops!"
   end do
   x1_d = x1
   call ieee_get_flag(ieee_usual, fpe_flag)
   lerr = lerr .or. any(fpe_flag)  
   if (lerr) write(*,*) "BFGS did not converge!"
   call ieee_set_status(original_fpe_status)

   return 

   contains

      function gradf(f, N, x1, dx) result(grad)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! Purpose:  Estimates the gradient of a function using a central difference
         !! approximation
         !! Inputs:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   N    :  number of variables N
         !!   x1   :  x value array
         !!   dx   :  step size to use when calculating derivatives
         !! Returns
         !!   grad :  N sized array containing estimated gradient of f at x1
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)  :: N
         class(lambda_obj),      intent(in)  :: f
         real(QP), dimension(:), intent(in)  :: x1
         real(QP),               intent(in)  :: dx
         ! Result
         real(QP), dimension(N)              :: grad
         ! Internals
         integer(I4B) :: i, j
         real(DP), dimension(N) :: xp, xm
         real(DP) :: fp_d, fm_d, dx_d
         real(QP) :: fp, fm

         dx_d = dx
         do i = 1, N
            do j = 1, N
               if (j == i) then
                  xp(j) = x1(j) + dx_d
                  xm(j) = x1(j) - dx_d
               else
                  xp(j) = x1(j)
                  xm(j) = x1(j)
               end if
            end do
            fp_d = f%eval(xp)
            fm_d = f%eval(xm)
            fp = real(fp_d, kind=QP)
            fm = real(fm_d, kind=QP)
            grad(i) = (fp - fm) / (2 * dx)
         end do
         return 
      end function gradf

      function minimize1D(f, x0, S, N, eps, lerr) result(astar)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This program find the minimum of a function of N variables in a single direction
         !! S using in sequence:
         !!    1.  A Bracketing method
         !!    2.  The golden section method
         !!    3.  A quadratic polynomial fit
         !! Inputs
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   N    :  Number of variables of function f
         !!   eps  :  Accuracy of 1 - dimensional minimization at each step
         !! Output
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! Returns
         !!   astar      :  Final minimum along direction S
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)  :: N
         class(lambda_obj),      intent(in)  :: f
         real(QP), dimension(:), intent(in)  :: x0, S
         real(QP),               intent(in)  :: eps
         logical,                intent(out) :: lerr
         ! Result
         real(QP)                            :: astar
         ! Internals
         integer(I4B) :: num = 0
         real(QP), parameter :: step = 0.7_QP     !! Bracketing method step size   
         real(QP), parameter :: gam = 1.2_QP      !! Bracketing method expansion parameter   
         real(QP), parameter :: greduce = 0.2_QP  !! Golden section method reduction factor   
         real(QP), parameter :: greduce2 = 0.1_QP ! Secondary golden section method reduction factor   
         real(QP) :: alo, ahi                     !! High and low values for 1 - D minimization routines   
         real(QP), parameter :: a0 = epsilon(1.0_QP)       !! Initial guess of alpha   
        
         alo = a0
         call bracket(f, x0, S, N, gam, step, alo, ahi, lerr)
         if (lerr) then
            write(*,*) "BFGS bracketing step failed!"
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if
         call golden(f, x0, S, N, greduce, alo, ahi, lerr)
         if (lerr) then
            write(*,*) "BFGS golden section step failed!"
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if
         call quadfit(f, x0, S, N, eps, alo, ahi, lerr)
         if (lerr) then
            write(*,*) "BFGS quadfit failed!"
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if 
         ! Quadratic fit method won't converge, so finish off with another golden section   
         call golden(f, x0, S, N, greduce2, alo, ahi, lerr)
         if (.not. lerr) astar = (alo + ahi) / 2.0_QP
         return 
      end function minimize1D

      function n2one(f, x0, S, N, a) result(fnew)
         implicit none
         ! Arguments
         integer(I4B),           intent(in) :: N
         class(lambda_obj),      intent(in) :: f
         real(QP), dimension(:), intent(in) :: x0, S
         real(QP),               intent(in) :: a
         ! Return
         real(QP) :: fnew
         ! Internals
         real(DP), dimension(N) :: xnew
         integer(I4B) :: i
         
         xnew(:) = real(x0(:) + a * S(:), kind=DP)
         fnew = real(f%eval(xnew(:)), kind=QP)
         return 
      end function n2one

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      subroutine bracket(f, x0, S, N, gam, step, lo, hi, lerr)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This subroutine brackets the minimum.  It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   gam  :  expansion parameter
         !!   step :  step size
         !!   lo   :  initial guess of lo bracket value
         !! The outputs include
         !!   lo   :  lo bracket
         !!   hi   :  hi bracket
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(in)    :: f
         real(QP), dimension(:), intent(in)    :: x0, S
         real(QP),               intent(in)    :: gam, step
         real(QP),               intent(inout) :: lo
         real(QP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals
         real(QP) :: a0, a1, a2, a3, da
         real(QP) :: f0, f1, f2, fcon
         integer(I4B) :: i
         integer(I4B), parameter :: MAXLOOP = 1000 ! maximum number of loops before method is determined to have failed   
         real(QP), parameter :: eps = epsilon(lo) ! small number precision to test floating point equality   
         real(QP), parameter :: dela = 12.4423_QP ! arbitrary number to test if function is constant   

         ! set up initial bracket points   
         lerr = .false.
         a0 =  lo
         da = step
         a1 = a0 + da
         a2 = a0 + 2 * da
         f0 = n2one(f, x0, S, N, a0)
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
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
            else if ((f0 < f1) .and. (f1 < f2)) then ! function appears to increase   
               da = da * gam
               a3 = a0 - da
               a2 = a1
               a1 = a0
               a0 = a3
               f2 = f1
               f0 = n2one(f, x0, S, N, a0)
            else if ((f0 < f1) .and. (f1 > f2)) then ! more than one minimum present, so in this case we arbitrarily choose the RHS min   
               da = da * gam
               a3 = a2 + da
               a0 = a1
               a1 = a2
               a2 = a3
               f0 = f1
               f1 = f2
               f2 = n2one(f, x0, S, N, a2)
            else if ((f0 > f1) .and. (abs(f2 - f1) <= eps)) then ! RHS equal   
               da = da * gam
               a3 = a2 + da
               a2 = a3
               f2 = n2one(f, x0, S, N, a2)
            else if ((abs(f0 - f1) < eps) .and. (f2 > f1)) then ! LHS equal   
               da = da * gam
               a3 = a0 - da
               a0 = a3
               f0 = n2one(f, x0, S, N, a0)
            else  ! all values equal stops if there is no minimum or takes RHS min if it exists   
               ! test if function itself is constant   
               fcon = n2one(f, x0, S, N, a2 + dela) !add by an arbitrary number to see if constant  
               if (abs(f2 - fcon) < eps) then   
                  lerr = .true.
                  return ! function is constant   
               end if
               a3 = a0 + 0.5_QP * (a1 - a0)
               a0 = a1
               a1 = a2
               a2 = a3
               f0 = f1
               f1 = f2     
               a3 = a0 + 0.5_QP * (a1 - a0)
               a1 = a2
               a2 = a3
               f1 = f2
               f2 = n2one(f, x0, S, N, a2)
            end if
         end do
         lerr = .true.
         return ! no minimum found   
      end subroutine bracket

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      subroutine golden(f, x0, S, N, eps, lo, hi, lerr) 
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses the golden section method to reduce the starting interval lo, hi by some amount sigma.  
         !! It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   gam  :  expansion parameter
         !!   eps  :  reduction interval in range (0 < sigma < 1) such that:
         !!             hi(new) - lo(new) = eps * (hi(old) - lo(old))
         !!   lo   :  initial guess of lo bracket value
         !! The outputs include
         !!   lo   :  lo bracket
         !!   hi   :  hi bracket
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(in)    :: f
         real(QP), dimension(:), intent(in)    :: x0, S
         real(QP),               intent(in)    :: eps
         real(QP),               intent(inout) :: lo
         real(QP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals 
         real(QP), parameter :: tau = 0.5_QP * (sqrt(5.0_QP) - 1.0_QP)  ! Golden section constant   
         integer(I4B), parameter :: MAXLOOP = 40 ! maximum number of loops before method is determined to have failed (unlikely, but could occur if no minimum exists between lo and hi)   
         real(QP) :: i0 ! Initial interval value   
         real(QP) :: a1, a2
         real(QP) :: f1, f2
         integer(I4B) :: i, j

         i0 =  hi - lo
         a1 =  hi - tau * i0
         a2 =  lo + tau * i0
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
         lerr = .false.
         do i = 1, MAXLOOP 
            if (abs((hi - lo) / i0) <= eps) return ! interval reduced to input amount   
            if (f2 > f1) then
               hi = a2
               a2 = a1
               f2 = f1
               a1 = hi - tau * (hi - lo)
               f1 = n2one(f, x0, S, N, a1)
            else 
               lo = a1
               a1 = a2
               f2 = f1
               a2 = hi - (1.0_QP - tau) * (hi - lo)
               f2 = n2one(f, x0, S, N, a2)
            end if
         end do
         lerr = .true.
         return ! search took too many iterations - no minimum found   
      end subroutine golden

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

      subroutine quadfit(f, x0, S, N, eps, lo, hi, lerr) 
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses a quadratic polynomial fit to locate the minimum of a function
         !! to some accuracy eps.  It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   lo    :  low bracket value
         !!   hi    :  high bracket value
         !!   eps   :  desired accuracy of final minimum location
         !! The outputs include
         !!   lo   :  final minimum location
         !!   hi   :  final minimum location
         !! Notes: Uses the ieee_exceptions intrinsic module to allow for graceful failure due to floating point exceptions, which won't terminate the run.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(in)    :: f
         real(QP), dimension(:), intent(in)    :: x0, S
         real(QP),               intent(in)    :: eps
         real(QP),               intent(inout) :: lo
         real(QP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals 
         integer(I4B), parameter :: MAXLOOP = 20 ! maximum number of loops before method is determined to have failed.   
         real(QP) :: a1, a2, a3, astar   ! three points for the polynomial fit and polynomial minimum   
         real(QP) :: f1, f2, f3, fstar   ! three function values for the polynomial and polynomial minimum   
         real(QP), dimension(3) :: row_1, row_2, row_3, rhs, soln        ! matrix for 3 equation solver (gaussian elimination)   
         real(QP), dimension(3,3) :: lhs
         real(QP) :: d1, d2, d3, aold, denom, errval
         integer(I4B) :: i

         lerr = .false.
         ! Get initial a1, a2, a3 values   
         a1 =  lo
         a2 =  lo + 0.5_QP * (hi - lo)
         a3 =  hi
         aold = a1
         astar = a2
         f1 = n2one(f, x0, S, N, a1)
         f2 = n2one(f, x0, S, N, a2)
         f3 = n2one(f, x0, S, N, a3)
         do i = 1, MAXLOOP 
            ! check to see if convergence is reached and exit   
            errval = abs((astar - aold) / astar)
            call ieee_get_flag(ieee_usual, fpe_flag)
            if (any(fpe_flag)) then
               lerr = .true.
               exit
            end if
            if (errval < eps) then
               lo = astar
               hi = astar
               exit
            end if
            ! Set up system for gaussian elimination equation solver   
            row_1 = [1.0_QP, a1, a1**2]
            row_2 = [1.0_QP, a2, a2**2]
            row_3 = [1.0_QP, a3, a3**2]
            rhs = [f1, f2, f3]
            lhs(1, :) = row_1
            lhs(2, :) = row_2
            lhs(3, :) = row_3
            ! Solve system of equations   
            soln(:) = util_solve_linear_system(lhs, rhs, 3, lerr)
            call ieee_set_flag(ieee_all, .false.) ! Set all flags back to quiet
            if (lerr) exit
            aold = astar
            if (soln(2) == soln(3)) then ! Handles the case where they are both 0. 0/0 is an unhandled exception
               astar = 0.5_QP
            else
               astar =  -soln(2) / (2 * soln(3))
            end if
            call ieee_get_flag(ieee_usual, fpe_flag)
            if (any(fpe_flag)) then
               lerr = .true.
               exit
            end if
            fstar = n2one(f, x0, S, N, astar)
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
         if (lerr) return
         lo = a1
         hi = a3
         return 
      end subroutine quadfit

end function util_minimize_bfgs