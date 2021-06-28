submodule(swiftest_data_structures) s_user_get_token
contains
   module function user_get_token(buffer, ifirst, ilast, ierr) result(token)
      !! author: David A. Minton
      !!
      !! Retrieves a character token from an input string. Here a token is defined as any set of contiguous non-blank characters not 
      !! beginning with or containing "!". If "!" is present, any remaining part of the buffer including the "!" is ignored
      !!
      !! Adapted from David E. Kaufmann's Swifter routine user_get_token.f90
      use swiftest, except_this_one => user_get_token
      implicit none

      ! Arguments
      character(len=*), intent(in)     :: buffer         !! Input string buffer
      integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
      integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
      integer(I4B), intent(out)        :: ierr           !! Error code
      character(len=:),allocatable     :: token          !! Returned token stringn

      ! Internals
      integer(I4B) :: i,ilength

      ilength = len(buffer)

      if (ifirst > ilength) then
         ilast = ifirst
         ierr = -1 !! Bad input
         token = ''
         return
      end if
      do i = ifirst, ilength
         if (buffer(i:i) /= ' ') exit
      end do
      if ((i > ilength) .or. (buffer(i:i) == '!')) then
         ifirst = i
         ilast = i
         ierr = -2 !! No valid token
         token = ''
         return
      end if
      ifirst = i
      do i = ifirst, ilength
         if ((buffer(i:i) == ' ') .or. (buffer(i:i) == '!')) exit
      end do
      ilast = i - 1
      ierr = 0

      token = buffer(ifirst:ilast)

      return

   end function user_get_token
end submodule s_user_get_token
