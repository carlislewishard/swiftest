pure subroutine util_crossproduct(ar1, ar2, ans)
   !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Calculates cross product of two arrays. Stores intermediate values in quad precision to improve accuracy when the two input
   !! arrays are near alignment
   use swiftest_globals
   use module_interfaces, EXCEPT_THIS_ONE => util_crossproduct
   implicit none
   ! Arguments
   real(DP),dimension(:),intent(in)  :: ar1,ar2
   real(DP),dimension(:),intent(out) :: ans
   ! Internals
   real(QP), dimension(3) :: qar1, qar2, qans 

   qar1(:) = real(ar1(:), kind=QP)
   qar2(:) = real(ar2(:), kind=QP)
   qans(1) = qar1(2) * qar2(3) - qar1(3) * qar2(2)
   qans(2) = qar1(3) * qar2(1) - qar1(1) * qar2(3)
   qans(3) = qar1(1) * qar2(2) - qar1(2) * qar2(1)
   ans(:) = real(qans(:), kind=DP)
   
  return

end subroutine util_crossproduct