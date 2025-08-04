      module fluxg
      contains
       function
      parameter(lmax=800, mmax=800)
      double precision :: aa(5,5,5),vova,vv(0:2,1:8,0:lmax,0:mmax)
!       intent(in) :: vv,step
!       intent(out) :: aa,bb
!       contains 
!         module pressg
!         double precision :: ek,gamma,vv(),uu(),beta,gamma 
! !kinetic energy density
!          ek(:,:) = 0.5d0 * gamma * (vv(2,step,:,:)**2+vv(3,step,:,:)**2)/vv(1,step,:,:)
! !internal energy density
!          eint(:,:) = vv(5,step,:,:)-ek(:,:) 
! !note that eint is really eint*rho but you never need the true value of e anyway
!          ptot(:,:) =  0.5*(beta*(gamma-1)*eint(:,:))
!          save ptot
!         end module pressg
!       end contains
 
 
       
!The flux functions for ideal MHD:

! Flux in x direction, ff. third dimension commented out for now.
! drho/dtx = rho * dvx  /dx
       aa(1,:,:) = vova*vv(2,step,:,:)
!dvx/dtx = rho*vx^2 + 1/gamma Ptot /dx
       aa(2,:,:) = vova*((vv(2,step,:,:)**2/vv(1,step,:,:)) + ptot(:,:)) 

!dvy/dtx = rho* vxvy - bxby
       aa(3,:,:) = vova*(vv(2,step,:,:)*vv(3,step,:,:)/vv(1,step,:,:))
!dvz/dtx = rho*vxvz - bxbz
!      ff(4,:,:) = vova*vv(2,step,:,:)*vv(4,step,:,:)/vv(1,step,:,:))
!       aa(4,:,:) = 0d0
!dE/dtx  = (E+Ptot)vx 
       aa(5,:,:)=  vova*((vv(5,step,:,:)+ptot(:,:))*vv(2,step,:,:)/vv(1,step,:,:))


!Flux in the Y direction, gg. z terms commented out for now 

!drho/dty =  rho*vy
       bb(1,:,:) = vova*vv(3,step,:,:)
!dvx/dty  = rho(v1v2) 
       bb(2,:,:) = vova*(vv(2,step,:,:)*vv(3,step,:,:)/vv(1,step,:,:))
!dvy/dty  = rho(v2v2) - b2b2 + ptot
       bb(3,:,:) = vova*((vv(3,step,:,:)**2)/vv(1,step,:,:) + ptot(:,:))
!dvz/dty  = rho(v3v2) - b2b3
!       bb(4,:,:) = vv(4,:,:))*vv(3,:,:)/vv(1,:,:) - vv(6,:,:)*vv(8,:,:)
!      bb(4,:,:) = 0
!de/dty   = (e+ptot)v2 - b2(b.v)
       bb(5,:,:) = vova*((vv,step(5,:,:)+ptot(:,:))*vv,step(3,:,:)/vv,step(1,:,:))
       
       save, aa, bb
!Flux in z direction, cc. left blank for now but easily inserted

      end module fluxg

 