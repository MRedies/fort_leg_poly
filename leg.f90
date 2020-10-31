program main
   use m_npy
   real    :: x(1000)
   real    :: step 
   integer :: n
   real, allocatable :: P(:,:)

   step = 2.0 / (size(x) - 1)

   x(1) = -1.0
   do n = 2,size(x)
      x(n) = x(n-1) + step
   enddo 

   call legendre_poly(10, x, P)

   call save_npy("leg.npy", P)

contains
   subroutine legendre_poly(n_max, x, P)
      implicit none
      integer, intent(in)              :: n_max
      real, intent(in)                 :: x(:)
      real, intent(inout), allocatable :: P(:,:)

      integer :: n

      if(allocated(P)) then 
         if(size(P,1) /= size(x) .or. size(P,2) /= n_max+1) then 
            deallocate(P)
         endif 
      endif 
      
      if(.not. allocated(P)) then
         allocate(P(size(x), 0:n_max), source=0.0)
      endif

      if (n_max >= 0) then
         P(:, 0) = 1.0
      endif

      if (n_max >= 1) then
         P(:, 1) = x
      endif

      do n = 1, n_max - 1
         P(:, n + 1) = ((2*n + 1)*x*P(:, n) - n*P(:, n - 1))/(n + 1.0)
      enddo
   end subroutine legendre_poly
end program main
