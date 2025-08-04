     Implicit none
      Integer&
      nt,lt,mt,l,m,n,nprint,lmax,mmax,lp,mp,lm,mm,i,p,&
      nr,ncentre,noff,time,mdiff,msum,m1,m2,x,ionise
     
      Parameter(lmax=700, mmax=700)
      Character filnam(14)*9, code*1,fdamp,experiment*7
      Double Precision                                              &
      uu(1:8,0:2,0:lmax,0:mmax),ff(1:8,0:lmax,0:mmax),                  &
      gg(1:8,0:lmax,0:mmax),tt(1:8,0:lmax,0:mmax),                  &
      bdotv(0:lmax,0:mmax),pb(0:lmax,0:mmax),                             &
      eint(0:lmax,0:mmax),aa(1:8,0:lmax,0:mmax) , sigma(1:8,0:lmax,0:mmax),    &
      bb(1:8,0:lmax,0:mmax),ttg(1:5,0:lmax,0:mmax),  gasinit(1:5,0:lmax,0:mmax), plasinit(1:8,0:lmax,0:mmax),   &
      vv(1:5,0:2,0:lmax,0:mmax),vova, vvt(1:5,0:lmax,0:mmax),uut(1:8,0:lmax,0:mmax),   &
      etot(0:lmax,0:mmax),ptot(0:lmax,0:mmax),mu,prssr,vb,             &
      offset,decay,wamp,decx,damp,dect,kx,xterm,fr,zterm,pi2           &
      ,spbomb(1:201,1:199),ek(0:lmax,0:mmax),emin,vmag,vreln(2),ener,vrel(2),       &
      gamma,beta,eobo2,eb(0:lmax,0:mmax),et(0:lmax,0:mmax),aistart     &
      ,frac,dpodg,bxo,byo,eo,pelletrate,boundmin,ioni,spbombamp,fracai
     
      namelist/datin/lt,mt,nt,nprint,mu,nr,noff,wamp,decx,kx,fr,experiment,spbombamp,aistart        &
      ,dect,m1,m2,damp,code,beta,frac,dpodg,bxo,byo,pelletrate,boundmin,fdamp,ionise,emin,fracai
      common/vars/ff,gg,uu,etot,ptot,vb,bdotv,eobo2,beta,aa,bb,vv,pb,bxo,byo,eo,vova,gasinit,sigma,vvt,plasinit,uut
      common/handy/filnam,mu,lt,mt,nt,nprint,noff,wamp                 &
      ,decx,damp,m1,m2,dect,kx,fr,pi2,code,spbomb,experiment                     &
      ,gamma,frac,dpodg,tt,pelletrate,boundmin,emin,fracai,ionise,aistart
