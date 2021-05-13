module lambda_function
   !! Defines a class that can enable objects that behave like lambda functions.
   use swiftest_globals
   implicit none

   type, public :: lambda_obj 
      private
      procedure(lambda0), pointer, nopass :: lambdaptr => null()
   contains
      generic   :: init => lambda_init_0
      procedure :: eval => lambda_eval_0
      procedure :: lambda_init_0
   end type

   abstract interface
      function lambda0(x)
         ! Template for a 0 argument function
         import DP
         real(DP), intent(in) :: x
         real(DP) :: lambda
      end function
   end interface

   contains
      subroutine lambda_init_0(self, lambda)
         implicit none
         ! Arguments
         class(lambda_obj), intent(out) :: self
         procedure(lambda0)             :: lambda
   
         self%lambdaptr => lambda
      end subroutine lambda_init_0
   
      function lambda_eval_0(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj), intent(in) :: self
         real(DP),          intent(in) :: x
         ! Result
         real(DP)                      :: y
   
         if (associated(self%lambdaptr)) then
            y = self%lambdaptr(x)
         else
            error stop "Initialize the object (call init) before computing values (call exec)!"
         end if
      end function lambda_eval_0

end module lambda_function

