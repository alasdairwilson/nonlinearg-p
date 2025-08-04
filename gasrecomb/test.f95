program test
double precision a(1:3),b(1:3),c(1:3)

   a(1) = (1)
   a(2) = (2)
   a(3) = (3)

   b(1) = (1)
   b(2) = (2)
   b(3) = (3)

   c(1) = (0)
   c(2) = (0)
   c(3) = (0)

where (a*b.gt.4) c = 1

print*, c(1:3)

end program
