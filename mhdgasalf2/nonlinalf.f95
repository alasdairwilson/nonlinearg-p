      Module declarations
      contains
      Include 'nonlinalf.inc'
      end module
      
      Program nonling
      use declarations 
      

! This program is a non-linear finite difference solver for use in numerical MHD-gas interactions. The program sets initial conditions in initialise, has a finite difference solver in solve (which can be switched out to various schemes) and a function for calculation the flux functions in flux (in all 3 dimensions) for MHD. Flux functions for gas are in Fluxg. Data is written out by function ritout. The program is exited in wrapup. Initial plan is to write it in 3-d but with 3-d functionality disabled, to see how well it runs before tackling harder problems.
!---------------------------------------------------------------
!our variables in MHD are stored in tensor uu(*,-,-,-,-):
! 1 - density
! 2 - vx
! 3 - vy
! 4 - vz
! 5 - energy
! 6 - bx
! 7 - by
! 8 - bz


!for the solve current time is indexed  second in uu(-,*,-,-,-)
!0 - current, n
!1 - predictor half step
!2 - next, n+1

!our variables for the GAS are stored in tensor vv(*,-,-,-,-)
! 1 - Density
! 2 - vx
! 3 - vy
! 4 - vz
! 5 - energy

!the time and space indexing are the same as for the MHD

!For 2d and 3d the dimensions are next in uu (-,-,*,*,*) for x,y,z respectively
!---------------------------------------------------------------
      call initialise
      
      call solve

      call wrapup

      end
!---------------------------------------------------------------

      Subroutine initialise
       
      Include 'nonlinalf.inc'
!   set handy constants  
! 
       pi2 = 2.d0*dacos(-1.d0)

       gamma=5.d0/3.d0

       filnam(1)='plde .out'
       filnam(2)='plvx .out'
       filnam(3)='plvy .out'
       filnam(4)='plvz .out'
       filnam(5)='plen .out'
       filnam(6)='plbx .out'
       filnam(7)='plby .out'
       filnam(8)='plbz .out'
       
       filnam(9)='gade .out'
       filnam(10)='gavx .out'
       filnam(11)='gavy .out'
       filnam(12)='gavz .out'
       filnam(13)='gaen .out'
!First load in variables from datin file.
       open(10,file='nonlinalf.in', status='old')

       read(10,datin)
       write(6,datin)
           
      do i=1,13
       write(filnam(i)(5:5),'(a1)') code
       open(20+i,file=filnam(i), status='unknown')
       rewind(20+i)
      end do

       write(6,*) 'Input complete: ','lt=',lt,'mt=',mt

!In this subroutine the matrices of initial values are created. 

!check lt,mt don't exceed array dimensions. I wonder if there is a better way to have the array size set in program rather than in the datin file to avoid this.
      if(lt.ge.lmax) then
       lt=lmax-1
       print*,'Warning: lt reset to ',lt
      end if
       if(mt.ge.mmax) then
       mt=mmax-1
      print*,'Warning: mt reset to ',mt
      end if
      do i = 0,lt+1
       dummym(i)=i
      end do
      ! E0 = et + eb

! E0/B0^2 = 1/2 * [ beta / gamma-1 + 1 ]
!set our two sound speed scaling constants, one is the ratio of magnetic energy to total energy which is directly related to beta, it is our scaling factor for the plasma energy
        eobo2 = 0.5d0 * ( beta/(gamma-1) + 1 )
        print*,eobo2
!the other is our scaling factor of the gas, it is again a function of beta and is the ratio of gas sound speed to alfven speed.
!The derivation of this expression comes from the assumption that the gas abd plasma are at the same temperature.
!vo/va = (gamma*beta)^0.5
        vova = (gamma*beta)**0.5d0


!next the setup of initial conditions for the plasma
   


!everything set to 0
        
         uu(:,:,:,:) = 0.0d0
        
!initial conditions set next
!density, sensible to set rho = 1, ie the density is equal to the normalized density everywhere, density gradients can be included.

        uu(1,:,:,:) = 1.0d0
!velocities set to 0 initially.
        uu(5,:,:,:)= 1.0d0
!Bx = can set any magnetic field, Obviously make sure initial magnetic field meets the requirements of divB = 0. The code is NOT divB conserving so it will slowly drift from a true solution to the magnetic field.

        uu(6,:,:,:) = bxo
        uu(7,:,:,:) = byo
!internal energy per unit mass (function of temperature)
        eint(:,:)   = 1.0d0
        
        !we set mean conditions for the sacrificial region. this is ugly better to read in initial conditions from a file and set or to calculate the mean of each variable (mean or equilibrium position).
! it is also entirely possible to make it a 1d array with no spatial dependence if the mean conditions are the same throughout the computational domain. Saves computational time at the expense for a more complex extrap function.

!first the plasma
      plasinit(1,:,:)   = 1.d0
      plasinit(2:4,:,:) =  0.d0
      plasinit(5,:,:)   = 1.0d0
      plasinit(6,:,:)   = 1.0d0
      plasinit(7:8,:,:) = 0.0d0
!and the gas
      gasinit(1,:,:)    = 1.d0
      gasinit(2:4,:,:)  = 0.0d0
      gasinit(5,:,:)    =  1.0d0
!Very important to ensure that initial values are consistent with the chosen plasma parameters (b0 v0 E0 Beta etc.)
!       beta is read in but our other normalization constant E0/B0^2 is a function of beta.


!Next we can add initial conditions that deviate from our equilibrium:

!--------------------------------------------------------------------------------------------------------
!!This code segment sets up a gaussian of density in the middle of the plasma from file spbomb.in
      open(101,file='spbomb.in', status='old')
       do l=1,201
       read(101,101)(spbomb(l,m),m=1,199)
       end do 
      open(501,file='spbomb.out',status='unknown')
      do l = 1,201
      write(501,101)(spbomb(l,m),m=1,199)
      end do
 
!        do i = 0,2     
!          uu(1,i,225:375,225:375) = (1.0d0+spbombamp*spbomb(25:175,25:175))
!        end do
!          uu(5,:,:,:)=uu(1,:,:,:)
!--------------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------------
! sigma, the damping factor. it may not error but you will get anomoulous results.
!remember that values returned by the code for any cell outside of the roi (ie with sigma>0.0) will be UNPHYSICAL. The damping layer is purely to prevent the boundary from influencing the roi.   
       print*, 'asig = ',asigma, 'bsig = ', bsigma
       do l=0,200
         ldum(l) = asigma*dexp(-(bsigma*dble(l)))
         ldum(lt-l) = ldum(l)
       end do
       do m = 0,200

        mdum(m) = asigma*dexp(-(bsigma*dble(m)))         
        mdum(mt-m) =  mdum(m)
       end do

       do l = 0,lt
        do m = 1,lt
         do i =1,8
         sigma(i,l,m) = ldum(l)+mdum(m)
        end do
       end do
       end do        
       print*,sigma(1,50,300),sigma(1,100,300),sigma(1,200,300),sigma(1,300,300)
!--------------------------------------------------------------------------------------------------------



!MUST DO ENERGY LAST, IT DEPENDS ON THE OTHER INTIAL VALUES
!setting the energy density is subtle E = Ek + Eb + Et, sum of the contributions from kinetic energy and magnetic pressure are easy but the last term, the gas pressure where p = (gamma-1)(density*energyperunitmass) requires assumptions. The simplest is that the temperature is constant and therefore the pressure is proportional to the density:


!the energy was initially set to 1.0d0 but if you have added a disturbance to the equilibrium at t=0 then the energy must adjust accordingly.
!with initial velocity field being 0 everywhere we assume that E0 = Eb+Ethermal then it is easy to express E as E0+Ek
!        ek(:,:) = 0.5d0*(uu(1,0,:,:)*((uu(2,0,:,:)+uu(3,0,:,:))/uu(1,0,:,:)))**2     
!        uu(1,0,:,:) = uu(1,0,:,:)
!        et = (beta/(2d0*(gamma-1)))*uu(1,0,:,:)*eint(:,:)
! 
!        eb = 0.5d0*(uu(6,0,:,:)**2+uu(7,0,:,:)**2)
! 
!        etot(:,:) = (1/eobo2)*(ek + et + eb)
!        
!        
!        uu(5,0,:,:) = etot(:,:)
!calculation of energy at equilibrium for fixed boundaries asssumes vx=vy=vz=0
!        eo = (1/eobo2)*(beta/(2d0*(gamma-1))*1.d0*1.d0+0.5d0*(bxo**2+byo**2))
! 
!        Print*, eo
!And for the gas must set ratio of densities from input files for coupling and vova *(v over va) which is vgassound/valfven which is equal to sqrt(gamma*Beta)


      print*, 'Beta = ',beta,'vova=',vova

      vv(:,:,:,:) = 0.0d0
        
!initial conditions set next
!density, sensible to set rho = 1, ie the density is equal to the normalized density everywhere
       vv(1,:,:,:) = 1.d0
!       vv(1,0,225:375,225:375) = 1.0d0 + spbombamp*spbomb(25:175,25:175)   
       
!there is no magnetic field to set for the gas

!next comes the energy calculation, similar to the one for the plasma with the absence of the amgnetic field  contribution to total energy density.
!since E0 = P0 for a cell with v = 0 then E = rho
      vv(5,:,:,:) = vv(1,:,:,:)


      print*,'vvvv',vv(1,0,200,3)

101    format(300(f8.6))        

      

        
      call ritout(0)


      end

!--------------------------------------------------------------

      Subroutine Solve
       
      Include 'nonlinalf.inc'
!This subroutine contains the finite difference solver. The richtmyer two-step scheme in 3 dimensions is used, this can be swapped out for any scheme that also only requires the flux functions, call the flux functions.

      do n = 1,nt

!looping over timesteps

       if(n.gt.200) then

	do m = 0,(mt+1)
	 do l = 0,(lt+1)
	    

         
! Kinetic energy present in the relative velocity.
	  vrel(1) = uu(3,l,m,0)-vv(3,l,m,0)
	  
	  vmag = dsqrt(vrel(1)**2)
	  ener = 0.5d0*vmag**2

	  
          
! calculate if energy present > threshold
	if(ionise.eq.1) then	
          if(ener.gt.emin) then
	   
  	 
           
         
 	   ioni = ((ener-emin)*frac)

! 	   
!            if(m.eq.0.and.l.eq.0) then
!              vreln(1) = dsqrt(vmag**2 - (2*ioni))
!              vreln(2) = (vreln(1)/vmag)
!              vel = vel*vreln(2)
!            end if         
!            
	   
	    uu(1,l,m,0) =  uu(1,l,m,0) + ioni
            if(l.eq.95.and.m.eq.95)print*,'ionising', ioni, vmag
        vv(1,l,m,0) =  vv(1,l,m,0) - ioni
 	  
           
! Add recombination possibly not needed, timescales are different.

! need to sutbract our used energy from the plasma velocity.

	   
 	   vreln(1) = dsqrt(vmag**2 - (2*ioni))
 

 	   vreln(2) = (vreln(1)/vmag)

 	   uu(2,0,l,m) =uu(2,0,l,m)*vreln(2)
	   uu(3,0,l,m) =uu(3,0,l,m)*vreln(2)
 	   vv(2,0,l,m) =vv(2,0,l,m)*vreln(2)
 	   vv(3,0,l,m) =vv(3,0,l,m)*vreln(2)
	   
! constituency check, is the magnitude of vrel the same now as vreln(1).

           end if
           end if
           end do
          end do
       end if


!calculate flux functions for all l,m at t = n
!not sure if necessary but set temp array to be current step, will reuse for corrector step, can also write two flux functions accessing different time indices of uu will see.
!First the plasma
     call driver(n)
      uut(:,:,:) =  uu(:,0,:,:)
      

      call flux
!Next the lax predictor step itself we use t =  1 for the predicted half-step

   

        uu(1:8,1,1:lt,1:mt) = 0.25d0*((uu(1:8,0,2:lt+1,1:mt) + uu(1:8,0,0:lt-1,1:mt) + &
     uu(1:8,0,1:lt,2:mt+1) + uu(1:8,0,1:lt,0:mt-1))) - &
     0.25d0*mu*( ff(1:8,2:lt+1,1:mt) + gg(1:8,1:lt,2:mt+1) - ff(1:8,0:lt-1,1:mt) - gg(1:8,1:lt,0:mt-1))
!              
      

!call the flux function  on the half timestep 
   
       uut(:,:,:) =  uu(:,1,:,:)
!      print*,tt(1:8,95,95)
     
       call extrap

       call flux


!next the leapfrog corrector step to update uu to the full timestep  

!loop the variables and spatial dimensions


         uu(1:8,2,1:lt,1:mt)  = uu(1:8,0,1:lt,1:mt) - 0.5d0*mu*(ff(1:8,2:lt+1,1:mt) - &
       ff(1:8,0:lt-1,1:mt) + gg(1:8,1:lt,2:mt+1) - gg(1:8,1:lt,0:mt-1))

! next we copy the new timestep to the current in order for us to start over again. 

   
       uut(:,:,:) =  uu(:,2,:,:)
       call extrap
 
  
!        print*, uu(1,2,100,100), uu(6,2,100,1), uu(7,2,100,100)
       uu(:,0,:,:) = uut(:,:,:)



!----------------------------------------------------------------------------------------------
!Next, almost identical setup for the gas.     

      vvt(:,:,:) =  vv(:,0,:,:)


      call fluxg
!Next the lax predictor step itself we use t =  1 for the predicted half-step

!DEBUG      print*, aa(1,90,90),aa(2,90,90),aa(3,90,90),aa(5,90,90)

!loop the variables and spatial dimensions


        vv(1:5,1,1:lt,1:mt) = 0.25d0*((vv(1:5,0,2:lt+1,1:mt) + vv(1:5,0,0:lt-1,1:mt) + &
       vv(1:5,0,1:lt,2:mt+1) + vv(1:5,0,1:lt,0:mt-1))) - &
       0.25d0*mu*( aa(1:5,2:lt+1,1:mt) + bb(1:5,1:lt,2:mt+1) - aa(1:5,0:lt-1,1:mt) - bb(1:5,1:lt,0:mt-1))
             


   

!call the flux function  on the half timestep 

      vvt(:,:,:) =  vv(:,1,:,:)

      call extrapg

      call fluxg
 

!next the leapfrog corrector step to update uu to the full timestep  

!loop the  spatial dimensions

         vv(1:5,2,1:lt,1:mt)  = vv(1:5,0,1:lt,1:mt) - 0.5d0*mu*(aa(1:5,2:lt+1,1:mt) - &
       aa(1:5,0:lt-1,1:mt) + bb(1:5,1:lt,2:mt+1) - bb(1:5,1:lt,0:mt-1))


!-------------------------------------------------------------------------------------
!Using the Macormack scheme with [FxBy, BxFy] on the [predictor, corrector]
! 
!       do n = 1,nt
! !      call driver(n)
! !solving for the gas
!   
!       ttg(:,:,:) = vv(:,0,:,:)
! 
!       call fluxg(ttg)
! 
! !the predictor step is forward in x, backward in y
! 
!       do i = 1,5
! !loop the variables and spatial dimensions
! 	do l = 1,lt
! 	  do m = 1, mt
! 
! !useful values lplus lminus etc.
! 	  lp = l + 1
! 	  lm = l - 1
!           mp = m + 1
!           mm = m - 1
! 
!         vv(i,1,l,m) = 0.25d0*((vv(i,0,lp,m) + vv(i,0,lm,m) + &
!        vv(i,0,l,mp) + vv(i,0,l,mm))) - &
!         mu*( aa(i,lp,m) + bb(i,l,m) - aa(i,l,m) - bb(i,l,mm))
!              
!            end do
!          end do   
!        end do
! 
! 
! !call the flux function  on the predictor 
! 
!       ttg(:,:,:) =  vv(:,1,:,:)
! 
!       call extrapg
!       call fluxg(ttg)
!      
! 
! !next the corrector step has the form backward in x, forward in y. note it uses half the timestep unlike the predictor 
!       do i = 1,5
! !loop the variables and spatial dimensions
! 	do l = 1,lt
! 	  do m = 1, mt
! 	  lp = l + 1
!           lm = l - 1
!           mp = m + 1
!           mm = m - 1
! 
! 	    vv(i,2,l,m)  = 0.5d0*(vv(i,0,l,m)+vv(i,1,l,m)) - & 
! 	    0.5d0*mu*(aa(i,l,m) - aa(i,lm,m) + bb(i,l,mp) - bb(i,l,m))
!           end do
!         end do
!        end do
! -------------------------------------------------------------------------------

! next we copy the new timestep to the current in order for us to start over again. 
       vvt(:,:,:) =  vv(:,2,:,:)
       call extrapg
!        print*, uu(1,2,100,100), uu(6,2,100,1), uu(7,2,100,100)
       vv(:,0,:,:) = vvt(:,:,:)      

!-------------------------------------------------------------------------------
! !Next we do the momentum coupling. 
!to disable coupling set frac to 0
!the subtelety here is in (dpodg/vova)
!for the gas in x.
! 
        vv(2,0,:,:) = vv(2,0,:,:)+frac*((dpodg/vova)*uu(2,0,:,:)-vv(2,0,:,:))

!y      
        vv(3,0,:,:) = vv(3,0,:,:)+frac*((dpodg/vova)*uu(3,0,:,:)-vv(3,0,:,:))
!for the plasma
!x

        uu(2,0,:,:) = uu(2,0,:,:)+frac*((vova/dpodg)*vv(2,0,:,:)-uu(2,0,:,:))

!y
        uu(3,0,:,:) = uu(3,0,:,:)+frac*((vova/dpodg)*vv(3,0,:,:)-uu(3,0,:,:))
!output every nprint steps
! 
      if(mod(n,nprint).eq.0) call ritout(n)
!       print*,'debug',uu(1,0,300,300)
      end do
      end
!-------------------------------------------------------------

      Subroutine Flux
       
      Include 'nonlinalf.inc'

!This subroutines calculates (up to 3) flux functions based on the suuplied values of uu. For the plasma.
!FF = xprint*, size(uut,1)
!GG = y
!HH = z if needed
      

!first calculate energy pressure and v.b for use in the flux functions.
       call press
 
       bdotv(:,:)= uut(2,:,:)*uut(6,:,:) + uut(3,:,:)*uut(7,:,:)
       
!The flux functions for ideal MHD:

! Flux in x direction, ff. third dimension commented out for now.
! drho/dtx = rho * dvx  /dx
       ff(1,:,:) = uut(2,:,:)
!dvx/dtx = rho*vx^2 - bx^2 + Ptot /dx 
       ff(2,:,:) = (uut(2,:,:)**2/uut(1,:,:)) - uut(6,:,:)**2 &
        + Ptot(:,:)
!dvy/dtx = rho* vxvy - bxby
       ff(3,:,:) =  uut(2,:,:)*uut(3,:,:)/uut(1,:,:) - uut(6,:,:)*uut(7,:,:)
!dvz/dtx = rho*vxvz - bxbz
!      ff(4,:,:) = uut(2,:,:)*uut(4,:,:)/uut(1,:,:)) - uut(6,:,:)*uut(8,:,:)
       ff(4,:,:) = 0d0
!dE/dtx  = (E+Ptot)vx - bx*(b.v)      
       ff(5,:,:)=(1d0/eobo2)*((uut(5,:,:)+ptot(:,:))*uut(2,:,:)/uut(1,:,:) &
       - uut(6,:,:)*bdotv(:,:))
!dbx/dtx = 0
       ff(6,:,:) = 0
!dby/dtx = vxby - vybx
       ff(7,:,:) =(uut(2,:,:)*uut(7,:,:)- uut(3,:,:)*uut(6,:,:))/uut(1,:,:)
!dbz/dtx = -(vybz - vzby)
!      ff(8,:,:) = -(uut(3,:,:)*uut(8,:,:) - t(4,:,:)*uut(7,:,:))/uut(1,:,:)
       ff(8,:,:)  = 0
       

!       print*, ff(2,0,100),ff(2,1,100),ff(2,2,100),ff(2,3,100)
!Flux in the Y direction, gg. z terms commented out for now 

!drho/dty =  rho*vy
       gg(1,:,:) = uut(3,:,:)
!dvx/dty  = rho(v1v2) - b1b2
       gg(2,:,:) = uut(2,:,:)*uut(3,:,:)/uut(1,:,:) - uut(6,:,:)*uut(7,:,:)
!dvy/dty  = rho(v2v2) - b2b2 + ptot
       gg(3,:,:) = (uut(3,:,:)**2)/uut(1,:,:) - uut(7,:,:)**2 + ptot(:,:)
!dvz/dty  = rho(v3v2) - b2b3
!       gg(4,:,:) = uut(4,:,:))*uut(3,:,:)/uut(1,:,:) - uut(6,:,:)*uut(8,:,:)
       gg(4,:,:) = 0
!de/dty   = (e+ptot)v2 - b2(b.v)
       gg(5,:,:) =(1/eobo2)*((uut(5,:,:)+ptot(:,:))*uut(3,:,:)/uut(1,:,:)  &
        - uut(7,:,:)*bdotv(:,:))
!dbx/dty  = -(v1b2-v2b1)
       gg(6,:,:) = -(uut(2,:,:)*uut(7,:,:) - uut(3,:,:)*uut(6,:,:))/uut(1,:,:)
!dby/dty  = 0
       gg(7,:,:) = 0
!dbz/dtz  = v2b3-v3b2
!      gg(8,:,:) = (uut(3,:,:)*uut(8,:,:) - uut(4,:,:)*uut(7,:,:))/uut(1,:,:)
       gg(8,:,:)  = 0

!Flux in z direction, hh. left blank for now but easily inserted

      end

!----------------------------------------------------------
      Subroutine Fluxg
       

      Include 'nonlinalf.inc'

!This subroutines calculates (up to 3) flux functions based on the suuplied values of uu.
!AA = x
!BB = y
!CC = z if needed
      

!first calculate energy pressure and v.b for use in the flux functions.
       call pressg
 
       
!The flux functions for ideal MHD:

! Flux in x direction, ff. third dimension commented out for now.
! drho/dtx = rho * dvx  /dx
       aa(1,:,:) = vova*vvt(2,:,:)
!dvx/dtx = rho*vx^2 + 1/gamma Ptot /dx
       aa(2,:,:) = vova*((vvt(2,:,:)**2/vvt(1,:,:))  + Ptot(:,:)) 

!dvy/dtx = rho* vxvy - bxby
       aa(3,:,:) =  vova*vvt(2,:,:)*vvt(3,:,:)/vvt(1,:,:)
!dvz/dtx = rho*vxvz - bxbz
!      ff(4,:,:) = vova*vvt(2,:,:)*vvt(4,:,:)/vvt(1,:,:))
       aa(4,:,:) = 0d0
!dE/dtx  = (E+Ptot)vx 
       aa(5,:,:)=vova*((vvt(5,:,:)+ptot(:,:))*vvt(2,:,:)/vvt(1,:,:))


!Flux in the Y direction, gg. z terms commented out for now \

!drho/dty =  rho*vy
       bb(1,:,:) = vova*vvt(3,:,:)
!dvx/dty  = rho(v1v2) 
       bb(2,:,:) = vova*vvt(2,:,:)*vvt(3,:,:)/vvt(1,:,:)
!dvy/dty  = rho(v2v2) - b2b2 + ptot
       bb(3,:,:) = vova*(vvt(3,:,:)**2)/vvt(1,:,:) + ptot(:,:)
!dvz/dty  = rho(v3v2) - b2b3
!       bb(4,:,:) = vvt(4,:,:))*vvt(3,:,:)/vvt(1,:,:) - vvt(6,:,:)*vvt(8,:,:)
       bb(4,:,:) = 0
!de/dty   = (e+ptot)v2 - b2(b.v)
       bb(5,:,:) =vova*((vvt(5,:,:)+ptot(:,:))*vvt(3,:,:)/vvt(1,:,:))
       

!Flux in z direction, cc. left blank for now but easily inserted

      end
!-------------------------------------------------------------

      subroutine press
      Include 'nonlinalf.inc'

!this subroutine calculates the pressure of the plasma via adiabatic law.

!magnetic energy density 
      pb(:,:) = 0.5d0*(uut(6,:,:)**2+uut(7,:,:)**2)
!kinetic energy density
      ek(:,:) = 0.5d0 * uut(1,:,:) * ((uut(2,:,:)/uut(1,:,:))**2+  &
      (uut(3,:,:)/uut(1,:,:))**2)
!internal energy density
      eint(:,:) = uut(5,:,:)-0.5*pb(:,:)-ek(:,:) 
!note that eint is really eint*rho but you never need the true value of e anyway
      ptot(:,:) =  0.5*(beta*(gamma-1)*eint(:,:)+2.d0*pb(:,:))
! this is representing the equation as 
! it is of interest how this can be adjusted, P can be a tensor with cross terms, here it is just listed as a scalar. 
! the gas part of the pressure comes from P = nRT/v, eint is really proportional to T and since density is known and V is held constant pressure can be worked out using p = (gamma-1) * (rho*eint)   eint is proportional to RT. 
      end
!----------------------------------------------------------
      subroutine pressg
      Include 'nonlinalf.inc'

!pressure function for the gas.


!kinetic energy density
      ek(:,:) = 0.5d0 * gamma * (vvt(2,:,:)**2+vvt(3,:,:)**2)/vvt(1,:,:)
!internal energy density
      eint(:,:) = vvt(5,:,:)-ek(:,:) 
!note that eint is really eint*rho but you never need the true value of e anyway
      ptot(:,:) =  eint(:,:)
 
      end
!-------------------------------------

      Subroutine extrap
      include 'nonlinalf.inc'

!WITH SACRIFICIAL REGION, FIXED BOUNDARIES
     
       uut(1:8,1:lt,1:mt) = uut(1:8,1:lt,1:mt)-sigma(1:8,1:lt,1:mt)*(uut(1:8,1:lt,1:mt)-plasinit(1:8,1:lt,1:mt))  
       

!       integer nn
! ! !update the dummy points at the l=lt+1,m=0, m=mt+1 edges:
! ! ! use quadratic extrapolant
! 
!        do m=1,mt
! 
!         do i = 1,8
! !l=lt+1 case:
! 
!          tt(i,lt+1,m)=uut(i,lt-2,m)+3.d0*uut(i,lt,m)-3.d0*uut(i,lt-1,m)
! 
!         end do                              
!        end do
! !!l = 0 case
!       do m=1,mt
!         do i = 1,8
! !l=lt+1 case:
!          tt(i,0,m)=uut(i,3,m)+3.d0*uut(i,1,m)-3.d0*uut(i,2,m)
!         end do                              
!        end do
! 
! 
! !m=0
!        do l=1,lt
!         do i = 1,8
!      
!          tt(i,l,0)=uut(i,l,3)+3.d0*uut(i,l,1)-3.d0*uut(i,l,2)
!  
!         end do
!        end do
! 
! ! finally, the dummy points at the far end: m=mt+1
!        do l=1,lt
!         do i = 1,8
!        
!         tt(i,l,mt+1)=uut(i,l,mt-2)+3.d0*uut(i,l,mt)-3.d0*uut(i,l,mt-1)
! 
!         end do
!        end do
! !  
!end quadratic boundaries~

!For driving waves at the boundary, sensible to use the extrap subroutine for efficiency.




   

        end       
!---------------------------------------------------------
     Subroutine extrapg
      include 'nonlinalf.inc'





!A sacrifical damping layer with a quadratic extrapolant for the BC.

!first implement the damping layer
!       print*,vvt(1,200,3)


       vvt(1:5,1:lt,1:mt) = vvt(1:5,1:lt,1:mt)-sigma(1:5,1:lt,1:mt)*(vvt(1:5,1:lt,1:mt)-gasinit(1:5,1:lt,1:mt))  

!       print*,  vvt(1,200,3)
!this reduces the size of the fluctuations from initial conditions based on sigma. Sigma should be 0 in the region of interest and increase gradually through the sacrificial region.

!the end points are done as a quadratic extrapolant.

! ! update the dummy points at the l=lt+1,m=0, m=mt+1 edges:
! ! use quadratic extrapolant
! 
!        do m=1,mt
! 
!         do i = 1,5
! !l=lt+1 case:
! 
!          ttg(i,lt+1,m)=vvt(i,lt-2,m)+3.d0*vvt(i,lt,m)-3.d0*vvt(i,lt-1,m)
!          ttg(i,0,m)=vvt(i,3,m)+3.d0*vvt(i,1,m)-3.d0*vvt(i,2,m)
!         end do                              
!        end do
! !l = 0 case
! 
! ! m=0,mt+1 case:
!        do l=1,lt
!         do i = 1,5
! 
!          ttg(i,l,0)=vvt(i,l,3)+3.d0*vvt(i,l,1)-3.d0*vvt(i,l,2)
!          ttg(i,l,mt+1)=vvt(i,l,mt-2)+3.d0*vvt(i,l,mt)-3.d0*vvt(i,l,mt-1)
!         end do
! 
!        end do


! !A periodic Boundary, if I can't get this one to work then fucking hell who knows what to do next.
!         
! 
! 
! 
! 
!          if(n.lt.noff) then
!            do m = m1,m2
!            zterm=dsin(dble(n-1)*pi2*fr/dble(1+noff))
!            mdiff=m2-m1
!            msum=m2+m1
! 
! 
!             xterm = 1.d0
! !           xterm= dexp(-decx*(m-msum/2)**2)
! 
!             
!             ttg(2,1,m)=0.5d0*damp*1.d0*zterm*xterm 
!             
! !DEBUG LINE
!            if(m.eq.100) then
!             print*, 'driving', ttg(2,0,100)
!            end if
! !DEBUG END
!           end do
!           elseif(n.eq.noff) then
!            do m = 1,mt
!              ttg(2:4,lt+1,m) = 0.d0
!              ttg(2:4,0,m)    = 0.d0    
!            end do 
!           else
!            do m = 1,mt
!             ttg(:,lt+1,m) = ttg(:,1,m)
!             ttg(:,0,m)    = ttg(:,lt,m)    
!           end do 
!           end if
! 
!           do l=1,lt
!            ttg(:,l,mt+1)= ttg(:,l,1)
!            ttg(:,l,0)   = ttg(:,l,mt)
!           end do
! 





!A reflective boundary. Density and energy are constants (need not be updated) and velocities have their sign reversed at the boundary. DOESNT WORK

!      do m=1,mt
!       
!       ttg(2,lt+1,m) =  -1.d0*ttg(2,lt-1,m)
!       ttg(2:4,lt,m)   =  0.0d0
!       ttg(1,lt,m)     =  1.0d0
!       ttg(5,lt,m)     =  1.0d0
!       
!       ttg(2,0,m)    =  -1.d0*ttg(2,2,m)
!       ttg(2:4,1,m)    =  0.0d0
!       ttg(1,1,m)     =  1.0d0
!       ttg(5,1,m)     =  1.0d0
!      end do
! 
!      do l=1,lt
!       ttg(3,l,mt+1)   =  -1.d0*ttg(3,l,mt-1)
!       ttg(2:4,l,mt)   =   0.d0
!       ttg(1,l,mt)   =   1.d0
!       ttg(5,l,mt)   =   1.d0
!      
!       ttg(3,l,0)      =  -1.d0*ttg(3,l,2)
!       ttg(2:4,l,1)    =   0.d0
!       ttg(1,l,1)   =   1.d0
!       ttg(5,l,1)   =   1.d0 
!      end do

!set gradient at boundary to 0, ie dummy point equal to alst FD point

!     do m=1,mt
!     
!       ttg(:,lt+1,m) = ttg(:,lt,m)
!       ttg(:,0,m)    = ttg(:,1,m)
! 
!     end do
! 
! 
!     do l = 1,lt
!  
!       ttg(:,l,mt+1)=  ttg(:,l,mt)
!       ttg(:,l,0)=     ttg(:,l,1)
!  
!     end do
!----------------------------------------------------------------------------
!set gradient at boundary equal to internal gradient, a linear extrapolant	
!     do m=1,mt
!     
!       ttg(:,lt+1,m) = 2.d0*ttg(:,lt,m)-ttg(:,lt-1,m)
!       ttg(:,0,m)    = 2.d0*ttg(:,1,m)-ttg(:,2,m)
! 
!     end do
! 
! 
!     do l = 1,lt
!  
!       ttg(:,l,mt+1)=  2.d0*ttg(:,l,mt)-ttg(:,l,mt-1)
!       ttg(:,l,0)= 2.d0*ttg(:,l,1)-ttg(:,l,2)
!  
!     end do
!------------------------------------------------------------------------------------
!a  quadratic extrapolant


! update the dummy points at the l=lt+1,m=0, m=mt+1 edges:
! use quadratic extrapolant
! 
     
! 
!        
! !l=lt+1 cases
! 
!          where(abs(ttg(2:4,lt-2,1:mt)+3.d0*ttg(2:4,lt,1:mt)-3.d0*ttg(2:4,lt-1,1:mt)).lt.boundmin) 
!            ttg(2:4,lt+1,1:mt)=0.d0
!          elsewhere 
!            ttg(2:4,lt+1,1:mt)=ttg(2:4,lt-2,1:mt)+3.d0*ttg(2:4,lt,1:mt)-3.d0*ttg(2:4,lt-1,1:mt)
!          endwhere
! 
!          where(abs(ttg(1,lt-2,1:mt)+3.d0*ttg(1,lt,1:mt)-3.d0*ttg(1,lt-1,1:mt)-1.d0).lt.boundmin) 
!            ttg(1,lt+1,1:mt)=1.d0
!          elsewhere 
!            ttg(1,lt+1,1:mt)=ttg(1,lt-2,1:mt)+3.d0*ttg(1,lt,1:mt)-3.d0*ttg(1,lt-1,1:mt)
!          endwhere
! 
!          where(abs(ttg(5,lt-2,1:mt)+3.d0*ttg(5,lt,1:mt)-3.d0*ttg(5,lt-1,1:mt)-1.d0).lt.boundmin) 
!            ttg(5,lt+1,1:mt)=1.d0
!          elsewhere 
!            ttg(5,lt+1,1:mt)=ttg(5,lt-2,1:mt)+3.d0*ttg(5,lt,1:mt)-3.d0*ttg(5,lt-1,1:mt)
!          endwhere
! 
! !-----------------------------------------------------
! !l=0 cases:
!          where(abs(ttg(2:4,3,1:mt)+3.d0*ttg(2:4,1,1:mt)-3.d0*ttg(2:4,2,1:mt)).lt.boundmin)
!            ttg(2:4,0,1:mt) = 0.d0
!          elsewhere 
!            ttg(2:4,0,1:mt)=ttg(2:4,3,1:mt)+3.d0*ttg(2:4,1,1:mt)-3.d0*ttg(2:4,2,1:mt)
!          endwhere
!                  
!          where(abs(ttg(1,3,1:mt)+3.d0*ttg(1,1,1:mt)-3.d0*ttg(1,2,1:mt)-1.d0).lt.boundmin)
!            ttg(1,0,1:mt) = 1.d0
!          elsewhere 
!            ttg(1,0,1:mt)=ttg(1,3,1:mt)+3.d0*ttg(1,1,1:mt)-3.d0*ttg(1,2,1:mt)
!          endwhere        
! 
!          where(abs(ttg(5,3,1:mt)+3.d0*ttg(5,1,1:mt)-3.d0*ttg(5,2,1:mt)-1.d0).lt.boundmin)
!            ttg(5,0,1:mt) = 1.d0
!          elsewhere 
!            ttg(5,0,1:mt)=ttg(5,3,1:mt)+3.d0*ttg(5,1,1:mt)-3.d0*ttg(5,2,1:mt)
!          endwhere
! 
! !-------------------------------------------------------
! !m=0 cases:     
!          where(abs(ttg(2:4,1:lt,3)+3.d0*ttg(2:4,1:lt,1)-3.d0*ttg(2:4,1:lt,2)).lt.boundmin) 
!              ttg(2:4,1:lt,0) = 0.d0
!          elsewhere 
!              ttg(2:4,1:lt,0)=ttg(2:4,1:lt,3)+3.d0*ttg(2:4,1:lt,1)-3.d0*ttg(2:4,1:lt,2)
!          endwhere
! 
!          where(abs(ttg(1,1:lt,3)+3.d0*ttg(1,1:lt,1)-3.d0*ttg(1,1:lt,2)-1.d0).lt.boundmin) 
!              ttg(1,1:lt,0) = 1.d0
!          elsewhere 
!              ttg(1,1:lt,0)=ttg(1,1:lt,3)+3.d0*ttg(1,1:lt,1)-3.d0*ttg(1,1:lt,2)
!          endwhere
! 
!          where(abs(ttg(5,1:lt,3)+3.d0*ttg(5,1:lt,1)-3.d0*ttg(5,1:lt,2)-1.d0).lt.boundmin) 
!              ttg(5,1:lt,0) = 1.d0
!          elsewhere 
!              ttg(5,1:lt,0)=ttg(5,1:lt,3)+3.d0*ttg(5,1:lt,1)-3.d0*ttg(5,1:lt,2)
!          endwhere
! !-------------------------------------------------------
! ! m=mt+1 cases:          
!          where(abs(ttg(2:4,1:lt,3)+3.d0*ttg(2:4,1:lt,1)-3.d0*ttg(2:4,1:lt,2)).lt.boundmin) 
!            ttg(2:4,1:lt,0)=0.d0
!          elsewhere 
!            ttg(2:4,1:lt,mt+1)=ttg(2:4,1:lt,mt-2)+3.d0*ttg(2:4,1:lt,mt)-3.d0*ttg(2:4,1:lt,mt-1)
!          endwhere
! 
!          where(abs(ttg(1,1:lt,3)+3.d0*ttg(1,1:lt,1)-3.d0*ttg(1,1:lt,2)-1.d0).lt.boundmin) 
!            ttg(1,1:lt,0)=1.d0
!          elsewhere 
!            ttg(1,1:lt,mt+1)=ttg(1,1:lt,mt-2)+3.d0*ttg(1,1:lt,mt)-3.d0*ttg(1,1:lt,mt-1)
!          endwhere
! 
!          where(abs(ttg(5,1:lt,3)+3.d0*ttg(5,1:lt,1)-3.d0*ttg(5,1:lt,2)-1.d0).lt.boundmin) 
!            ttg(5,1:lt,0)=1.d0
!          elsewhere 
!            ttg(5,1:lt,mt+1)=ttg(5,1:lt,mt-2)+3.d0*ttg(5,1:lt,mt)-3.d0*ttg(5,1:lt,mt-1)
!          endwhere
! 
! ! 
!       
!-------------------------------------------------------------------------------------
! use quadratic extrapolant that tends to mean
! 
!        do m=1,mt
! 
!        
! !l=lt+1, l=0 cases:
! 
!          ttg(1,lt+1,m)  = (2.d0*(ttg(1,lt-2,m)+3.d0*ttg(1,lt,m)-3.d0*ttg(1,lt-1,m))+1.d0)/3.d0
!          ttg(5,lt+1,m)  = (2.d0*(ttg(5,lt-2,m)+3.d0*ttg(5,lt,m)-3.d0*ttg(5,lt-1,m))+1.d0)/3.d0
!          ttg(2:4,lt+1,m)= (2.d0*(ttg(2:4,lt-2,m)+3.d0*ttg(2:4,lt,m)-3.d0*ttg(2:4,lt-1,m)))/3.d0
! 
!          ttg(1,0,m)  = (2.d0*(ttg(1,3,m)+3.d0*ttg(1,1,m)-3.d0*ttg(1,2,m))+1.d0)/3.d0
!          ttg(5,0,m)  = (2.d0*(ttg(5,3,m)+3.d0*ttg(5,1,m)-3.d0*ttg(5,2,m))+1.d0)/3.d0
!          ttg(2:4,0,m)= (2.d0*(ttg(2:4,3,m)+3.d0*ttg(2:4,1,m)-3.d0*ttg(2:4,2,m)))/3.d0       
!                       
!        end do
! 
! 
!        do l=1,lt
! 
! !m=0, m=mt+1 cases:
!          ttg(:,l,0)=ttg(:,l,3)+3.d0*ttg(:,l,1)-3.d0*ttg(:,l,2)
!          
! 
!          ttg(1,l,0)=(2.d0*(ttg(1,l,3)+3.d0*ttg(1,l,1)-3.d0*ttg(1,l,2))+1.d0)/3.d0
!          ttg(5,l,0)=(2.d0*(ttg(5,l,3)+3.d0*ttg(5,l,1)-3.d0*ttg(5,l,2))+1.d0)/3.d0
!          ttg(2:4,l,0)=2.d0*(ttg(2:4,l,3)+3.d0*ttg(2:4,l,1)-3.d0*ttg(2:4,l,2))/3.d0
! 
!          ttg(1,l,mt+1)=(2.d0*(ttg(1,l,mt-2)+3.d0*ttg(1,l,mt)-3.d0*ttg(1,l,mt-1))+1.d0)/3.d0
!          ttg(5,l,mt+1)=(2.d0*(ttg(5,l,mt-2)+3.d0*ttg(5,l,mt)-3.d0*ttg(5,l,mt-1))+1.d0)/3.d0
!          ttg(2:4,l,mt+1)=2.d0*(ttg(2:4,l,mt-2)+3.d0*ttg(2:4,l,mt)-3.d0*ttg(2:4,l,mt-1))/3.d0 
!        end do
! !-----------------------------------------------------------------------------------------------------------
!quadratic extrapolant with last points excluded, or with corner adjacent points are linear
!       do m=2,mt-1
! 
!        
! !l=lt+1, l=0 cases:
! 
!          ttg(:,lt+1,m)=ttg(:,lt-2,m)+3.d0*ttg(:,lt,m)-3.d0*ttg(:,lt-1,m)
!          ttg(:,0,m)=ttg(:,3,m)+3.d0*ttg(:,1,m)-3.d0*ttg(:,2,m)                    
!        end do
! 
! 
!        do l=2,lt-1
! 
! !m=0, m=mt+1 cases:
!          ttg(:,l,0)=ttg(:,l,3)+3.d0*ttg(:,l,1)-3.d0*ttg(:,l,2)
!          ttg(:,l,mt+1)=ttg(:,l,mt-2)+3.d0*ttg(:,l,mt)-3.d0*ttg(:,l,mt-1)
! 
!        end do

!          ttg(:,lt+1,mt)=2.d0*ttg(:,lt,mt)-ttg(:,lt-1,mt)
!          ttg(:,lt,mt+1)=2.d0*ttg(:,lt,mt)-ttg(:,lt,mt-1)    
! 
!          ttg(:,lt+1,1)=2.d0*ttg(:,lt,1)-ttg(:,lt-1,1)
!          ttg(:,lt,0)  =2.d0*ttg(:,lt,1)-ttg(:,lt,2)   
! 
!          ttg(:,0,mt)=2.d0*ttg(:,1,mt)-ttg(:,2,mt)
!          ttg(:,1,mt+1)=2.d0*ttg(:,1,mt)-ttg(:,1,mt-1) 
! 
!          ttg(:,0,1)=2.d0*ttg(:,1,1)-ttg(:,2,1)
!          ttg(:,1,0)=2.d0*ttg(:,1,1)-ttg(:,1,2)   

! average of initial condition and current value
!

!         do m=1,mt
!         ttg(1,lt+1,m) =  (ttg(1,lt+1,m) + 1.0d0)/2
!         ttg(1,0,m)    =  (ttg(1,lt+1,m) + 1.0d0)/2
!  
!         ttg(5,lt+1,m) =  (ttg(5,lt+1,m) + 1.0d0)/2
!         ttg(5,0,m)    =  (ttg(5,0,m) + 1.0d0)/2
! 
!         ttg(2:3,lt+1,m) =  (ttg(2:3,lt+1,m) + 0.0d0)/2
!         ttg(2:3,0,m)    =  (ttg(2:3,0,m) + 0.0d0)/2
!         end do
! 
!         do l = 1,lt
! 
!         end do

        end       
!---------------------------------------------
        Subroutine driver(n)
        Include 'nonlinalf.inc'
!          if(n.lt.noff) then
!           tterm = dexp(-(n-noff/2d0)**2/1000d0)
!           uu(5,0,225:375,225:375) = uu(5,0,225:375,225:375)+0.001d0*spbomb(25:175,25:175)
!           vv(5,0,225:375,225:375) = vv(5,0,225:375,225:375)+0.001d0*spbomb(25:175,25:175)
!          endif
!        if(n.lt.noff) then
!             dummyl=dble(l)
!  
!         xterm(1:lt) = amp*dexp((-((dummym(1:lt))-decx)**2/35))
! !        yterm(1:20) = 1.d0*dexp((-((dummyl(1:20))-10d0)**2/9))
!         tterm = dsin(2*pi2*n/fr)
!        print*, 'balls'
! !  !           uut(1,222,201:399)= 0.2d0*xterm(201:399)*tterm+ uut(1,222,201:399)
! !           uut(1,221,201:399)= 0.5d0*xterm(201:399)*tterm+uut(1,221,201:399)
! !           uut(5,222,201:399)= 0.2d0*xterm(201:399)*tterm+uut(5,222,201:399)
! !           uut(5,221,201:399)= 0.5d0*xterm(201:399)*tterm+ uut(5,221,201:399)
! !         uut(1,218,201:399)= 0.2d0*xterm(201:399)*tterm+uut(1,218,201:399)
! !           uut(1,219,201:399)= 0.5d0*xterm(201:399)*tterm+uut(1,219,201:399)
! !           uut(5,218,201:399)= 0.2d0*xterm(201:399)*tterm+uut(5,218,201:399)
! !           uut(5,219,201:399)= 0.5d0*xterm(201:399)*tterm+uut(5,219,201:399)
! !           uut(1,220,201:399)= xterm(201:399)*tterm+uut(1,220,201:399)
! !           uut(5,220,201:399)= xterm(201:399)*tterm+uut(5,220,201:399)  
! 
!           uut(1,201:399,220)= xterm(201:399)*tterm+uut(1,220,201:399)
!           uut(5,201:399,220)= xterm(201:399)*tterm+uut(5,220,201:399)    
!   
! 
!         print*,uut(1,200,300),xterm(300),tterm
!         end if
!This driving term is for a gaussain shaped source for new gas for pellet ablation.  
!       if(n.lt.noff) then
!        ek(:,:) = 0.5d0 * gamma * ((vv(2,0,:,:)**2+vv(3,0,:,:)**2)/vv(1,0,:,:))
!             
!        vv(1,0,1:lt,1:mt) = vv(1,0,1:lt,1:mt) + pelletrate * spbomb(:,:)
! 
!        eint(:,:)= 1.d0
!       
!        et = vv(1,0,:,:)*eint(:,:)
!      
!        vv(5,0,:,:) = ek+et
!  
!       end if



      end


!----------------------------------------------------------
      Subroutine ritout(n)
       
      Include 'nonlinalf.inc'
!This routine writes out the current uu and vv values

      print*, 'Timestep =', n
      if(n.lt.5000) then
      

	write(21,*)'# timestep = ',n
	do l=200,400
	  write(21,100)((uu(1,0,l,m)-1.d0),',',m=200,400)
	end do
	write(21,*)
	write(21,*)
      
	write(22,*)'# timestep = ',n
	do l=200,400
	  write(22,100)(uu(2,0,l,m),',',m=200,400)
	end do
	write(22,*)
	write(22,*)

	write(23,*)'# timestep = ',n
	do l=200,400
	  write(23,100)(uu(3,0,l,m),',',m=200,400)
	end do
	write(23,*)
	write(23,*)

  !       write(24,*)'# timestep = ',n
  !       do l=200,400
  !       write(24,100)(uu(4,0,l,m),m=200,400)
  !       end do
  !       write(24,*)
  !       write(24,*)

	write(25,*)'# timestep = ',n
	do l=200,400
	  write(25,100)((uu(5,0,l,m)-1.d0),',',m=200,400)
	end do
	write(25,*)
	write(25,*)

	write(26,*)'# timestep = ',n
	do l=200,400
	  write(26,100)(uu(6,0,l,m),',',m=200,400)
	end do
	write(26,*)
	write(26,*)

	write(27,*)'# timestep = ',n
	do l=200,400
	  write(27,100)(uu(7,0,l,m),',',m=200,400)
	end do
	write(27,*)
	write(27,*)
    
	write(28,*)'# timestep = ',n
	do l=200,400
	  write(28,100)(uu(8,0,l,m),',',m=200,400)
	end do
	write(28,*)
	write(28,*)
    
	write(29,*)'# timestep = ',n
	do l=200,400
	  write(29,100)((vv(1,0,l,m)-1.d0),',',m=200,400)
	end do
	write(29,*)
	write(29,*)
    
	write(30,*)'# timestep = ',n
	do l=200,400
	  write(30,100)(vv(2,0,l,m),',',m=200,400)
	end do
	write(30,*)
	write(30,*)
    
	write(31,*)'# timestep = ',n
	do l=200,400
	  write(31,100)(vv(3,0,l,m),',',m=200,400)
	end do
	write(31,*)
	write(31,*)
    
	write(32,*)'# timestep = ',n
	do l=200,400
	  write(32,100)(vv(4,0,l,m),',',m=200,400)
	end do
	write(32,*)
	write(32,*)
    
	write(33,*)'# timestep = ',n
	do l=200,400
	  write(33,100)((vv(5,0,l,m)-1.d0),',',m=200,400)
	end do
	write(33,*)
	write(33,*)
      
      end if

	  
100   format(1000(f14.7,A1))


      end
!------------------------------------------------------------

      Subroutine wrapup
       
      Include 'nonlinalf.inc'
!this routine wraps up the code, it will likely do nothing, can be used to output current timestep to be used as input for a new run.


      end
      
!now for the functions      


!-------------------------------------------------------------    
