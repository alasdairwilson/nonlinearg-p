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
      gasinit(5,:,:)    =  gamma*1.0d0
!Very important to ensure that initial values are consistent with the chosen plasma parameters (b0 v0 E0 Beta etc.)
!       beta is read in but our other normalization constant E0/B0^2 is a function of beta.


!Next we can add initial conditions that deviate from our equilibrium:




!--------------------------------------------------------------------------------------------------------
! !this code reads in sigma, this is the damping factor. A file named in variable fdamp is loaded into sigma, the damping factor. !!THE SIZE OF FDAMP MUST MATCH THE GRIDSIZE EXACTLY it may not error but you will get anomoulous results.
! !remember that values returned by the code for any cell outside of the roi (ie with sigma>0.0) will be UNPHYSICAL. The damping layer is purely to prevent the boundary from influencing the roi.   
!       open(502,file='sigma200.in', status='old')
!        do l=1,lt
!        read(502,101)(sigma(1,l,m),m=1,mt)
!        end do 
!       open(503,file='sigma.out',status='unknown')
!       do l = 1,lt
!       write(503,101)(sigma(1,l,m),m=1,mt)
!       end do
! ! duplicate sigma for all variables not necessary since the switcfh to loop over variable instead of matrix multiplication
! !     sigma(2:8,l,m) = sigma(1,l,m)
!--------------------------------------------------------------------------------------------------------
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
     print*,'sigma', sigma(1,10,10), sigma(1,100,100), sigma(1,180,180), sigma(1,250,250)


 !      Print*, eo
!And for the gas must set ratio of densities from input files for coupling and vova *(v over va) which is vgassound/valfven which is equal to sqrt(gamma*Beta)


      print*, 'Beta = ',beta,'vova=',vova

      vv(:,:,:,:) = 0.0d0
        
!initial conditions set next
!density, sensible to set rho = 1, ie the density is equal to the normalized density everywhere
       vv(1,:,:,:) = 1.d0
       uu(1,:,:,:) = 1.d0 

print*, 'experiment chosen is  ',experiment
!-----------------------------------------------------------------------------------------------------------
     if(experiment.eq.1) then
       print*, 'Driver is spbomb in gas'
        !!This code segment sets up a the shape gaussian of density in the middle of the fluids from file spbomb.in
      open(101,file='spbomb.in', status='old')
       do l=1,201
       read(101,101)(spbomb(l,m),m=1,199)
       end do 
      open(501,file='spbomb.out',status='unknown')
      do l = 1,201
      write(501,101)(spbomb(l,m),m=1,199)
      end do

        
        vv(1,0,225:375,225:375) = 1.0d0 + spbombamp*spbomb(25:175,25:175)   
!-----------------------------------------------------------------------------------------------------------
elseif(experiment.eq.2) then
 print*, 'Driver is spbomb in plasma'
      open(101,file='spbomb.in', status='old')
       do l=1,201
       read(101,101)(spbomb(l,m),m=1,199)
       end do 
      open(501,file='spbomb.out',status='unknown')
      do l = 1,201
      write(501,101)(spbomb(l,m),m=1,199)
      end do
  
      uu(1,0,225:375,225:375) = 1.0d0 + spbombamp*spbomb(25:175,25:175)
elseif(experiment.eq.3) then
    print*, 'Driver is spbomb in both gas and plasma'
      open(101,file='spbomb.in', status='old')
       do l=1,201
       read(101,101)(spbomb(l,m),m=1,199)
       end do 
      open(501,file='spbomb.out',status='unknown')
      do l = 1,201
      write(501,101)(spbomb(l,m),m=1,199)
      end do


         vv(:,0,:,:)=gasinit
         uu(:,0,:,:)=plasinit

         uu(1,0,150:250,150:250) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,0,200:300,250:350) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,0,370:470,160:260) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,0,400:500,400:500) = (1.0d0-spbombamp*spbomb(50:150,50:150))
         uu(1,0,100:200,380:480) = (1.0d0-spbombamp*spbomb(50:150,50:150))

         vv(1,0,150:250,150:250) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,0,200:300,250:350) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,0,370:470,160:260) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,0,400:500,400:500) = (1.0d0+spbombamp*spbomb(50:150,50:150))
         vv(1,0,100:200,380:480) = (1.0d0+spbombamp*spbomb(50:150,50:150))
          
endif
       
!MUST DO ENERGY LAST, IT DEPENDS ON THE OTHER INTIAL VALUES
!setting the energy density is subtle E = Ek + Eb + Et, sum of the contributions from kinetic energy and magnetic pressure are easy but the last term, the gas pressure where p = (gamma-1)(density*energyperunitmass) requires assumptions. The simplest is that the temperature is constant and therefore the pressure is proportional to the density:


!the energy was initially set to 1.0d0 but if you have added a disturbance to the equilibrium at t=0 then the energy must adjust accordingly.
!with initial velocity field being 0 everywhere we assume that E0 = Eb+Ethermal then it is easy to express E as E0+Ek
       ek(:,:) = 0.5d0*(uu(1,0,:,:)*((uu(2,0,:,:)+uu(3,0,:,:))/uu(1,0,:,:)))**2     
       uu(1,0,:,:) = uu(1,0,:,:)
       et = (gamma*beta/(2d0*(gamma-1)))*uu(1,0,:,:)*eint(:,:)

       eb = 0.5d0*(uu(6,0,:,:)**2+uu(7,0,:,:)**2)

       etot(:,:) = (ek + et + eb)
       
!if onlythe density is changed then the energy gradients are equal to the density gradients:       
       uu(5,0,:,:) = gamma*uu(1,0,:,:)
!calculation of energy at equilibrium for fixed boundaries asssumes vx=vy=vz=0
!       eo = (1/eobo2)*(beta/(2d0*(gamma-1))*1.d0*1.d0+0.5d0*(bxo**2+byo**2))



!there is no magnetic field to set for the gas

!next comes the energy calculation, similar to the one for the plasma with the absence of the amgnetic field  contribution to total energy density.
!since E0 = P0 for a cell with v = 0 then E = rho
 
      vv(5,0,:,:) = gamma*vv(1,0,:,:)

      vv(:,1,:,:)=vv(:,0,:,:)
      vv(:,2,:,:)=vv(:,0,:,:)

      uu(:,1,:,:)=uu(:,0,:,:)
      uu(:,2,:,:)=uu(:,0,:,:)

     
!      print*,'vvvv',vv(1,0,200,3)

101    format(300(f8.6))        

      
       
        
      call ritout(0)


      end

!--------------------------------------------------------------

      Subroutine Solve
       
      Include 'nonlinalf.inc'
!This subroutine contains the finite difference solver. The richtmyer two-step scheme in 3 dimensions is used, this can be swapped out for any scheme that also only requires the flux functions, call the flux functions.
         
      do n = 1,nt

	  
      call driver(n)
     
!looping over timesteps

       if(n.gt.aistart) then
       if(ionise.eq.1) then
        if(n.eq.1)print*,'alfven ionization is enabled'
	do m = 0,(mt+1)
	 do l = 0,(lt+1)
          print*,'debug'
         
! Kinetic energy present in the relative velocity.
	  vrel(1) = uu(3,0,l,m)-vv(3,0,l,m)
         if(l.eq.350.and.m.eq.350) print*,vrel(1) ,uu(3,0,l,m),vv(3,0,l,m)
	  
	  vmag = dsqrt(vrel(1)**2)
	  ener = 0.5d0*vmag**2
          if(l.eq.350.and.m.eq.350) print*,ener ,emin
	  
          
! calculate if energy present > threshold
	
          if(ener.gt.emin) then
	   
  	 
           
         
 	   ioni = ((ener-emin)*fracai)

! 	   
!            if(m.eq.0.and.l.eq.0) then
!              vreln(1) = dsqrt(vmag**2 - (2*ioni))
!              vreln(2) = (vreln(1)/vmag)
!              vel = vel*vreln(2)
!            end if         
!            
	  
	    uu(1,0,l,m) =  uu(1,0,l,m) + ioni
            if(l.eq.355.and.m.eq.355)print*,'ionising', ioni, vmag
        vv(1,0,l,m) =  vv(1,0,l,m) - ioni
 	  

           
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
           end do
          end do
         end if
       end if



!calculate flux functions for all l,m at t = n
!not sure if necessary but set temp array to be current step, will reuse for corrector step, can also write two flux functions accessing different time indices of uu will see.
!First the plasma
     call driver(n)
      uut(:,:,:) =  uu(:,0,:,:)
      

      call flux
!Next the lax predictor step itself we use t =  1 for the predicted half-step
 
   

        uu(1:8,1,1:lt,1:mt) = 0.25d0*((uu(1:8,0,2:lp,1:mt) + uu(1:8,0,0:lt-1,1:mt) + &
     uu(1:8,0,1:lt,2:mt+1) + uu(1:8,0,1:lt,0:mt-1))) - &
     0.25d0*mu*( ff(1:8,2:lt+1,1:mt) + gg(1:8,1:lt,2:mt+1) - ff(1:8,0:lt-1,1:mt) - gg(1:8,1:lt,0:mt-1))
             
   

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
      eint(:,:) = uut(5,:,:)-pb(:,:)-ek(:,:) 
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
      ptot(:,:) =  0.5*(beta*(gamma-1)*eint(:,:))
 
      end
!-------------------------------------

      Subroutine extrap
      include 'nonlinalf.inc'

!WITH SACRIFICIAL REGION, FIXED BOUNDARIES
       
       do l = 1,mt
        do m = 1,lt
         do i = 1,8
          uut(i,l,m) = uut(i,l,m)-sigma(1,l,m)*(uut(i,l,m)-plasinit(i,l,m))
         end do
        end do
       end do  
       

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


     

!A sacrifical damping layer with fixed absolute boundaries. Signal should never reach the final points (0,lt+1) as the damping is sufficient.

!first implement the damping layer
!       print*,vvt(1,200,3)
!       print*,sigma(1,150,150),vvt(1,150,150)
       do l = 1,mt
        do m = 1,lt
         do i = 1,5
          vvt(i,l,m) = vvt(i,l,m)-sigma(1,l,m)*(vvt(i,l,m)-gasinit(i,l,m))
         end do
        end do
       end do  
!       print*,vvt(1,150,150)

!       print*,  vvt(1,200,3)
     

        end       
!---------------------------------------------
        Subroutine driver(n)
        Include 'nonlinalf.inc'

!while initial conditions are set in initialise a constant driving term can be added here.
! 
!       if(experiment.eq.'planewb') then
!       if(n.eq.1) print*,'planewl:plane wave driven from bottom boundary'
! !noff allows a turn off time for the driver.
!        if(n.lt.noff) then
! !m1 and m2 are the spatial extent for the driver,
!        do m=m1,m2
! 
! !* first the harmonic driver term..now we're driving perturbations against the magnetic
! !* field direction, that is, perturbations down the edge
! 
! !zterm is the sinousoidal term
!         zterm=dsin((dble(n-1)/500)*pi2*fr)
!         mdiff=m2-m1
!         msum=m2+m1
! 
! 
! ! next the spatial structure in the z (or m) direction, with gaussian
! ! envelope to smooth gradients
! !** we have a set of possibilities: either discrete sources, or a standing wave
! 
! !exponential envelope
!         xterm= dexp(-decx*(m-msum/2)**2)
! !plane wave
!         xterm = 1.d0
!    
! 
! !           tt(1,0,m)=1.d0+damp*1.d0*zterm*xterm
!           vv(2,0,198:200,m)=damp*1.d0*zterm*xterm
! 
!  
!        
!           uu(2,0,198:200,m)=damp*1.d0*zterm*xterm
! 
! 
! !           uu(1,1,0,m)=1.d0+damp*zterm*xterm*dexp(-dect*(n-noff))
! 
!         end do
! 
!         end if  
!        end if      



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


!      print*, vv(5,0,100,100)
      end


!----------------------------------------------------------
      Subroutine ritout(n)
       
      Include 'nonlinalf.inc'
          character(8)  :: date
    character(10) :: rtime
    character(5)  :: zone
    integer,dimension(8) :: values
    double precision :: timediff,nremaining
!This routine writes out the current uu and vv values

      print*, 'Timestep =', n

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

      write(21,*)'# timestep = ',n
      do l=0,lt+1
       write(21,100)((uu(1,0,l,m)-1.0d0),',',m=0,lt+1)
      end do
      write(21,*)
      write(21,*)

      write(22,*)'# timestep = ',n
      do l=0,lt+1
      write(22,100)(uu(2,0,l,m),',',m=0,lt+1)
      end do
       write(22,*)
       write(22,*)

      write(23,*)'# timestep = ',n
      do l=0,lt+1
      write(23,100)(uu(3,0,l,m),',',m=0,lt+1)
      end do
      write(23,*)
      write(23,*)

      write(24,*)'# timestep = ',n
      do l=0,lt+1
      write(24,100)(ptot(l,m),',',m=0,lt+1)
      end do
      write(24,*)
      write(24,*)

      write(25,*)'# timestep = ',n
      do l=0,lt+1
      write(25,100)((uu(5,0,l,m)-gamma*1.0d0),',',m=0,lt+1)
      end do
      write(25,*)
      write(25,*)

      write(26,*)'# timestep = ',n
      do l=0,lt+1
      write(26,100)(uu(6,0,l,m),',',m=0,lt+1)
      end do
      write(26,*)
      write(26,*)

      write(27,*)'# timestep = ',n
      do l=0,lt+1
      write(27,100)(uu(7,0,l,m),',',m=0,lt+1)
      end do
      write(27,*)
      write(27,*)
  
      write(28,*)'# timestep = ',n
      do l=0,lt+1
      write(28,100)(uu(8,0,l,m),',',m=0,lt+1)
      end do
      write(28,*)
      write(28,*)
  
      write(29,*)'# timestep = ',n
      do l=0,lt+1
      write(29,100)((vv(1,0,l,m)-1.d0),',',m=0,lt+1)
      end do
      write(29,*)
      write(29,*)
  
      write(30,*)'# timestep = ',n
      do l=0,lt+1
      write(30,100)(vv(2,0,l,m),',',m=0,lt+1)
      end do
      write(30,*)
      write(30,*)
  
      write(31,*)'# timestep = ',n
      do l=0,lt+1
      write(31,100)(vv(3,0,l,m),',',m=0,lt+1)
      end do
      write(31,*)
      write(31,*)
  
      write(32,*)'# timestep = ',n
      do l=0,lt+1
      write(32,100)(vv(4,0,l,m),',',m=0,lt+1)
      end do
      write(32,*)
      write(32,*)
  
      write(33,*)'# timestep = ',n
      do l=0,lt+1
      write(33,100)((vv(5,0,l,m)-gamma*1.0d0),',',m=0,lt+1)
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
