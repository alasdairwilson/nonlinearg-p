      Implicit none
      Integer
     & nt,lt,mt,l,m,n,nprint,lmax,mmax,lp,mp,lm,mm,i,p,
     & nr,ncentre,noff,time,mdiff,msum,m1,m2
     
      Parameter(lmax=250, mmax=250)
      Character filnam(10)*9, code*1
      Double Precision
     & uu(1:8,0:2,0:lmax,0:mmax),ff(1:8,0:lmax,0:mmax),
     & gg(1:8,0:lmax,0:mmax),tt(1:8,0:lmax,0:mmax),
     & bdotv(0:lmax,0:mmax),pb(0:lmax,0:mmax),u0,
     & eint(0:lmax,0:mmax),
     & etot(0:lmax,0:mmax),ptot(0:lmax,0:mmax),mu,prssr,vb,
     & offset,decay,wamp,decx,damp,dect,kx,xterm,fr,zterm,pi2
     & ,spbomb(1:202,1:202),kb,temp,ek(0:lmax,0:mmax),
     & gamma,beta,eobo2,eb(0:lmax,0:mmax),et(0:lmax,0:mmax)
     
      namelist/datin/lt,mt,nt,nprint,mu,nr,noff,wamp,decx,kx,fr
     & ,dect,m1,m2,damp,code,beta
      common/vars/ff,gg,uu,etot,ptot,vb,bdotv,eobo2,beta
      common/handy/filnam,mu,lt,mt,nt,u0,nprint,noff,wamp
     & ,decx,damp,m1,m2,dect,kx,fr,pi2,code,spbomb,kb,temp
     & ,gamma