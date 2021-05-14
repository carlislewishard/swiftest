module lambda_function
   !! Defines a class that can enable objects that behave like lambda functions.
   use swiftest_globals
   implicit none

   type, public :: lambda_obj 
      procedure(lambda0), pointer, nopass :: lambdaptr => null()

   contains
      generic   :: init => lambda_init_0
      procedure :: eval => lambda_eval_0
      procedure, nopass :: lambda_init_0
      final     :: lambda_destroy
   end type
   interface lambda_obj
      module procedure lambda_init_0
   end interface

   abstract interface
      function lambda0(x) result(y)
         ! Template for a 0 argument function
         import DP
         real(DP), dimension(:), intent(in) :: x
         real(DP)                           :: y
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
   
      function lambda_eval_0(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj),      intent(in) :: self
         real(DP), dimension(:), intent(in) :: x
         ! Result
         real(DP)                      :: y
   
         if (associated(self%lambdaptr)) then
            y = self%lambdaptr(x)
         else
            error stop "Lambda function was not initialized"
         end if
      end function lambda_eval_0

      subroutine lambda_destroy(self)
         implicit none
         type(lambda_obj) :: self
         if (associated(self%lambdaptr)) nullify(self%lambdaptr)
      end subroutine lambda_destroy

end module lambda_function

