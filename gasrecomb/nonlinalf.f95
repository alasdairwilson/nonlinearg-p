      Program nonling
       
      Include 'nonlinalf.inc'

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
      plasinit(5,:,:)   = gamma*1.0d0
      plasinit(6,:,:)   = 1.0d0
      plasinit(7:8,:,:) = 0.0d0
!and the gas
      gasinit(1,:,:)    = 1.d0
      gasinit(2:4,:,:)  = 0.0d0
      gasinit(5,:,:)    = gamma*1.0d0
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
 
       do i = 0,2     
         vv(:,i,:,:)=gasinit
         uu(:,i,:,:)=plasinit

         uu(1,i,150:250,150:250) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,i,200:300,250:350) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,i,370:470,160:260) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,i,400:500,400:500) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,i,100:200,380:480) = (1.0d0-spbombamp*spbomb(50:150,50:150))

         vv(1,i,150:250,150:250) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,i,200:300,250:350) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,i,370:470,160:260) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,i,400:500,400:500) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,i,100:200,380:480) = (1.0d0+spbombamp*spbomb(50:150,50:150))
       end do
         uu(5,:,:,:)=gamma*uu(1,:,:,:)
         vv(5,:,:,:)=gamma*vv(1,:,:,:)
!--------------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------------
! sigma, the damping factor. it may not error but you will get anomoulous results.
!remember that values returned by the code for any cell outside of the roi (ie with sigma>0.0) will be UNPHYSICAL. The damping layer is purely to prevent the boundary from influencing the roi.   
       print*, 'asig = ',asigma, 'bsig = ', bsigma
       do l=0,100
         ldum(l) = asigma*dexp(-(bsigma*dble(l)))
         ldum(lt-l) = ldum(l)
       end do
       do m = 0,100

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
       print*,sigma(1,50,300),sigma(1,100,300),sigma(1,200,300),sigma(1,300,300),sigma(1,300,500),sigma(1,300,550)
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
     
      uut(:,:,:) =  uu(:,0,:,:)
      call driver(n)

      call flux
!Next the lax predictor step itself we use t =  1 for the predicted half-step
      
      do i = 1,8
!loop the variables and spatial dimensions
	do l = 1,lt
	  do m = 1, mt
	lp = l + 1
	lm = l - 1
        mp = m + 1
        mm = m - 1

        uu(i,1,l,m) = 0.25d0*((uu(i,0,lp,m) + uu(i,0,lm,m) + &
     uu(i,0,l,mp) + uu(i,0,l,mm))) - &
     0.25d0*mu*( ff(i,lp,m) + gg(i,l,mp) - ff(i,lm,m) - gg(i,l,mm))
             
           end do
         end do   
       end do
   

!call the flux function  on the half timestep 
   
       uut(:,:,:) =  uu(:,1,:,:)
!      print*,tt(1:8,95,95)
     
       call extrap

       call flux


!next the leapfrog corrector step to update uu to the full timestep  
      do i = 1,8
!loop the variables and spatial dimensions
	do l = 1,lt
	  do m = 1, mt
	lp = l + 1
	lm = l - 1
        mp = m + 1
        mm = m - 1

         uu(i,2,l,m)  = uu(i,0,l,m) - 0.5d0*mu*(ff(i,lp,m) - &
       ff(i,lm,m) + gg(i,l,mp) - gg(i,l,mm))
          end do
        end do
       end do
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
      do i = 1,5
!loop the variables and spatial dimensions
	do l = 1,lt
	  do m = 1, mt
	lp = l + 1
	lm = l - 1
        mp = m + 1
        mm = m - 1

        vv(i,1,l,m) = 0.25d0*((vv(i,0,lp,m) + vv(i,0,lm,m) + &
       vv(i,0,l,mp) + vv(i,0,l,mm))) - &
       0.25d0*mu*( aa(i,lp,m) + bb(i,l,mp) - aa(i,lm,m) - bb(i,l,mm))
             
           end do
         end do   
       end do
   

!call the flux function  on the half timestep 

      vvt(:,:,:) =  vv(:,1,:,:)

      call extrapg

      call fluxg
 

!next the leapfrog corrector step to update uu to the full timestep  
      do i = 1,5
!loop the variables and spatial dimensions
	do l = 1,lt
	  do m = 1, mt
	lp = l + 1
	lm = l - 1
        mp = m + 1
        mm = m - 1

         vv(i,2,l,m)  = vv(i,0,l,m) - 0.5d0*mu*(aa(i,lp,m) - &
       aa(i,lm,m) + bb(i,l,mp) - bb(i,l,mm))
          end do
        end do
       end do




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
       
   

        end       
!---------------------------------------------------------
     Subroutine extrapg
      include 'nonlinalf.inc'

  





!first implement the damping layer
!       print*,vvt(1,200,3)
      

       vvt(1:5,1:lt,1:mt) = vvt(1:5,0:lt+1,0:mt+1)-sigma(1:5,0:lt+1,0:mt+1)*(vvt(1:5,0:lt+1,0:mt+1)-gasinit(1:5,0:lt+1,0:mt+1))  

!       print*,  vvt(1,200,3)
!this reduces the size of the fluctuations from initial conditions based on sigma. Sigma should be 0 in the region of interest and increase gradually through the sacrificial region.



        end       
!---------------------------------------------
        Subroutine driver(n)
        Include 'nonlinalf.inc'


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
          character(8)  :: date
    character(10) :: rtime
    character(5)  :: zone
    integer,dimension(8) :: values
    double precision :: timediff,nremaining
      if(n.eq.0) then
       call date_and_time(date,rtime,zone,values)
       t1=values
      end if
      call date_and_time(date,rtime,zone,values)
      t2=values
      timediff=60*60*t2(5)+60*t2(6)+t2(7)+0.001*t2(8)-60*60*t1(5)-60*t1(6)-t1(7)-0.001*t1(8)
      t1=t2
      nremaining=nt-n
      
  
      print*,'Estimated time remaining:', &
 floor((nremaining/20*timediff)/60), ' mins ', &
((nremaining/20*timediff)/60-floor((nremaining/20*timediff)/60))*60, &
 ' secs'
!This routine writes out the current uu and vv values 
      print*, 'Timestep =', n


      write(21,*)'# timestep = ',n
      do l=0,lt+1
       write(21,100)((uu(1,0,l,m)-1.d0),',',m=0,mt+1)
      end do
      write(21,*)
      write(21,*)

      write(22,*)'# timestep = ',n
      do l=0,lt+1
      write(22,100)(uu(2,0,l,m),',',m=0,mt+1)
      end do
       write(22,*)
       write(22,*)

      write(23,*)'# timestep = ',n
      do l=0,lt+1
      write(23,100)(uu(3,0,l,m),',',m=0,mt+1)
      end do
      write(23,*)
      write(23,*)

!       write(24,*)'# timestep = ',n
!       do l=0,lt+1
!       write(24,100)(uu(4,0,l,m),m=0,mt+1)
!       end do
!       write(24,*)
!       write(24,*)

      write(25,*)'# timestep = ',n
      do l=0,lt+1
      write(25,100)((uu(5,0,l,m)-1.d0),',',m=0,mt+1)
      end do
      write(25,*)
      write(25,*)

      write(26,*)'# timestep = ',n
      do l=0,lt+1
      write(26,100)(uu(6,0,l,m),',',m=0,mt+1)
      end do
      write(26,*)
      write(26,*)

      write(27,*)'# timestep = ',n
      do l=0,lt+1
      write(27,100)(uu(7,0,l,m),',',m=0,mt+1)
      end do
      write(27,*)
      write(27,*)
  
      write(28,*)'# timestep = ',n
      do l=0,lt+1
      write(28,100)(uu(8,0,l,m),',',m=0,mt+1)
      end do
      write(28,*)
      write(28,*)
  
      write(29,*)'# timestep = ',n
      do l=0,lt+1
      write(29,100)((vv(1,0,l,m)-1.d0),',',m=0,mt+1)
      end do
      write(29,*)
      write(29,*)
  
      write(30,*)'# timestep = ',n
      do l=0,lt+1
      write(30,100)(vv(2,0,l,m),',',m=0,mt+1)
      end do
      write(30,*)
      write(30,*)
  
      write(31,*)'# timestep = ',n
      do l=0,lt+1
      write(31,100)(vv(3,0,l,m),',',m=0,mt+1)
      end do
      write(31,*)
      write(31,*)
  
      write(32,*)'# timestep = ',n
      do l=0,lt+1
      write(32,100)(vv(4,0,l,m),',',m=0,mt+1)
      end do
      write(32,*)
      write(32,*)
  
      write(33,*)'# timestep = ',n
      do l=0,lt+1
      write(33,100)((vv(5,0,l,m)-1.d0),',',m=0,mt+1)
      end do
      write(33,*)
      write(33,*)
  

          
100   format(1000(f14.7,A1))


      end
!------------------------------------------------------------

      Subroutine wrapup
       
      Include 'nonlinalf.inc'
!this routine wraps up the code, it will likely do nothing, can be used to output current timestep to be used as input for a new run.


      end
      
!now for the functions      


!-------------------------------------------------------------    
