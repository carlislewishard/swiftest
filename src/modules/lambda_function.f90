module lambda_function
   !! Defines a class that can enable objects that behave like lambda functions.
   use swiftest_globals
   implicit none

   type, public :: lambda_obj 
      !! Base class for an lambda function object. This object takes no additional arguments other than the dependent variable x, an array of real numbers
      procedure(lambda0), pointer, nopass :: lambdaptr => null()
      real(DP) :: lastval
      real(DP),dimension(:), allocatable :: lastarg
   contains
      generic   :: init => lambda_init_0
      procedure :: eval => lambda_eval_0
      procedure, nopass :: lambda_init_0
      final     :: lambda_destroy
   end type

   type, public, extends(lambda_obj) :: lambda_obj_err
      !! Extended class for an lambda function object. This object takes allows for the return of a logical error flag during evaluation of the function.
      procedure(lambda0err), pointer, nopass :: lambdaptr_err => null()
      logical   :: lerr     
   contains
      generic   :: init => lambda_init_0_err
      procedure :: eval => lambda_eval_0_err
      procedure, nopass :: lambda_init_0_err
   end type
   interface lambda_obj
      module procedure lambda_init_0
      module procedure lambda_init_0_err
   end interface

   abstract interface
      function lambda0(x) result(y)
         ! Template for a 0 argument function
         import DP
         real(DP), dimension(:), intent(in) :: x
         real(DP)                           :: y
      end function
      function lambda0err(x, lerr) result(y)
         ! Template for a 0 argument function that returns an error value
         import DP
         real(DP), dimension(:), intent(in)  :: x
         logical,                intent(out) :: lerr
         real(DP)                            :: y
      end function
   end interface

   contains
      type(lambda_obj) function lambda_init_0(lambda)
         implicit none
         ! Arguments
         procedure(lambda0)             :: lambda
   
         lambda_init_0%lambdaptr => lambda
         return
      end function lambda_init_0

      type(lambda_obj_err) function lambda_init_0_err(lambda, lerr)
         implicit none
         ! Arguments
         procedure(lambda0err)  :: lambda
         logical, intent(in) :: lerr
         lambda_init_0_err%lambdaptr_err => lambda
         lambda_init_0_err%lerr = lerr
         return
      end function lambda_init_0_err
   
      function lambda_eval_0(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj),      intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         ! Result
         real(DP)                      :: y
   
         if (associated(self%lambdaptr)) then
            y = self%lambdaptr(x)
            self%lastval = y
            if (allocated(self%lastarg)) deallocate(self%lastarg)
            allocate(self%lastarg, source=x)
         else
            error stop "Lambda function was not initialized"
         end if
      end function lambda_eval_0

      function lambda_eval_0_err(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj_err),  intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         ! Result
         real(DP)                      :: y
   
         if (associated(self%lambdaptr_err)) then
            y = self%lambdaptr_err(x, self%lerr)
            self%lastval = y
            if (allocated(self%lastarg)) deallocate(self%lastarg)
            allocate(self%lastarg, source=x)
         else
            error stop "Lambda function was not initialized"
         end if
      end function lambda_eval_0_err

      subroutine lambda_destroy(self)
         implicit none
         type(lambda_obj) :: self
         if (associated(self%lambdaptr)) nullify(self%lambdaptr)
      end subroutine lambda_destroy

end module lambda_function

