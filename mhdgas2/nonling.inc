      Implicit none
      Integer&
      nt,lt,mt,l,m,n,nprint,lmax,mmax,lp,mp,lm,mm,i,p,&
      nr,ncentre,noff,time,mdiff,msum,m1,m2
     
      Parameter(lmax=250, mmax=250)
      Character filnam(14)*9, code*1
      Double Precision                                              &
      uu(1:8,0:2,0:lmax,0:mmax),ff(1:8,0:lmax,0:mmax),                  &
      gg(1:8,0:lmax,0:mmax),tt(1:8,0:lmax,0:mmax),                  &
      bdotv(0:lmax,0:mmax),pb(0:lmax,0:mmax),                             &
      eint(0:lmax,0:mmax),aa(1:5,0:lmax,0:mmax) ,                      &
      bb(1:5,0:lmax,0:mmax),ttg(1:5,0:lmax,0:mmax),                  &
      vv(1:5,0:2,0:lmax,0:mmax),vova,                                  &
      etot(0:lmax,0:mmax),ptot(0:lmax,0:mmax),mu,prssr,vb,             &
      offset,decay,wamp,decx,damp,dect,kx,xterm,fr,zterm,pi2           &
      ,spbomb(1:202,1:202),ek(0:lmax,0:mmax),                          &
      gamma,beta,eobo2,eb(0:lmax,0:mmax),et(0:lmax,0:mmax)        &
      ,frac,dpodg,bxo,byo,eo
     
      namelist/datin/lt,mt,nt,nprint,mu,nr,noff,wamp,decx,kx,fr        &
      ,dect,m1,m2,damp,code,beta,frac,dpodg,bxo,byo
      common/vars/ff,gg,uu,etot,ptot,vb,bdotv,eobo2,beta,aa,bb,vv,pb,bxo,byo,eo,vova
      common/handy/filnam,mu,lt,mt,nt,nprint,noff,wamp                 &
      ,decx,damp,m1,m2,dect,kx,fr,pi2,code,spbomb                      &
      ,gamma,frac,dpodg,ttg,tt