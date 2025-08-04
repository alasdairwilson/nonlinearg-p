       program tester
       integer lt, mt
       double precision spbomb(20,20)
       lt = 20
       mt = 20
       open(101,file='spbombsmall.in', status='old')
       rewind(101)
101    format(300(f8.6))        
       do l=1,lt
       read(101,101)(spbomb(l,m),m=1,mt)
       end do
100    format(300(f8.6))
       write(*,101)spbomb(10,:)
       end