 

Character filnam(14)*9
Double precision gasinit(1:5,0:lmax,0:mmax), plasinit(1:8,0:lmax,0:mmax),sigma(1:8,0:lmax,0:mmax) &
    ff(1:8,0:lmax,0:mmax),gg(1:8,0:lmax,0:mmax),&
    uu(1:8,0:2,0:lmax,0:mmax),vv(1:8,0:2,0:lmax,0:mmax),&
    vb,eobo2,beta,bdotv(0:lmax,0:mmax),pb(0:lmax,0:mmax) 


common/boundaries/gasinit,plasinit,sigma

common/vars/ff,gg,uu,vb,bdotv,eobo2,beta,aa,bb,vv,pb,bxo,byo,eo,vova

common/handy/filnam,mu,lt,mt,nt,nprint,noff,                &
      pi2,                    &
      ,gamma,frac,dpodg,emin,t1,t2
