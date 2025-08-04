       Program sigmacreate
 
       double precision ldum(0:600), mdum(0:600),asigma,bsigma,sigma(600,600)
       integer j,xx,yy
       
       xx = 600
       yy = 600
       asigma = 4.0d-1
       bsigma = 4.0d-2

       sigma(:,:) = 0.0d0

       print*,
       print*, 'asig = ',asigma, 'bsig = ', bsigma
       do i=0,200
         ldum(i) = asigma*dexp(-(bsigma*dble(i)))
         ldum(xx-i) = ldum(i)
       end do
       do j = 0,200
        mdum(j) = asigma*dexp(-(bsigma*dble(j)))
        mdum(yy-j) =  mdum(j)
       end do

       do i = 1,xx
        do j = 1,yy
         sigma(i,j) = sqrt((ldum(i))**2+(mdum(j))**2)
        end do
       end do        


 
       print*, sigma(100,100),sigma(101,100),sigma(102,100), sigma(250,250)
       



       
       open(501,file='sigma200.in',status='unknown')
      do l = 1,xx
      write(501,101)(sigma(l,m),m=1,yy)
      end do 
101   format(300(f8.6))   
       end


        