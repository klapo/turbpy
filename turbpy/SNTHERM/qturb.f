c***********************************************************************
c  QTURB returns wind functions qlat and qsen for the turbulent
c  transfer of latent and sensible heat. Streamlined version of qturb4.f 
c  created on 23 April,2001.
c***********************************************************************
      subroutine qturb(height,Tkair,Tsurface,Eair,Esurface,wsp,
     1wspo,qlat,qsen,dliqvol,cdryair,Rw,Ck,Csk,snowdepth,
     2ltype,z0,Cdn,rhoair,bpress,dlogT,dlogQ,dlogW)
c
c  Called from: MAIN(10),(16)
c
c  Arguments
c
c height(ld) : height above ground of measured met values
c Tkair : Air temperature [K]
c Tsurfaceo : Past surface temperature [K] when icall = 1
c           : Actual surface temperature [K] when icall = 2
c Eair : Vapor pressure of water in air [mb]
c Esurfaceo : Past vapor pressure at surface [mb] when icall = 1
c           : Actual vapor pressure at surface [mb] when icall = 2
c wsp :  Wind speed [m/s]
c wspo : Past wind speed [m/s]
c qlat : Wind function (coefficient) for latent heat exchange [W-m/kg]
c qlat0 : Same as above, except without windless term and not multipled 
c        by the latent heat of sublimation
c qsen : Wind function (coefficient) for sensible heat exchange 
c        [W/(m^2 K)]
c qsen0: Same as above, except without windless term 
c dliqvol : Volume fraction of liquid water
c cdryair : Specific heat of dry air at 0C (1005.0) [J/kg-k]
c Rw : gas constant for water vapor (461.296) [J/kg-K]
c Ck : Windless exchange coefficient for water vapor
c Csk : Windless exchange coefficient for heat
c snowdepth : snow depth [m]
c ltype(ld) :Layer type code  (1 = snow)
c z0 : Roughness length for turbulent transfer of momentum [m]
c Cdn : Drag coefficient at neutral stability
c rhoair : Density of air  [kg/m^3]
c bpress: Barometric pressure (mb)
c dlogT : ln(temperature height above surface/roughness length for heat)
c dlogQ : ln(humidity height above surface/roughness length for water vapor)
c dlogW : ln(wind helght above surface/roughness length for momentum)
c melt : =1 if in meltzone, otherwise = 0
c wstar : friction velocity  [m/s]
c icall : =1 for call to estimate qsen and qlat before heat balance
c         =2 for call after heat balance, to finalize transfer coefficients
c Rb : Bulk Richardson number
c
      double precision height(3),Tkair,Tsurface,Eair,Esurface,wsp
      double precision wspo,qlat,qsen,dliqvol
      double precision cdryair,Rw,Ck,Csk,snowdepth,z0,Cdn
      double precision rhoair,bpress,dlogT,dlogQ,dlogW
      integer ltype
c
c  Local
c 
c Hr: height above snowpack for relative humidity measurement
c      (Hr=height(3)-snowdepth).
c Ht: height above snowpack for temperature measurement
c      (Ht=height(1)-snowdepth).
c Hw: height above snowpack for wind speed measurement
c      (Hw=height(2)-snowdepth).
c Next 3 variables are the integrated stability corrections, all computed using
c   actual measurement heights
c Psim: Correction for momentum
c Psiq: Correction for water vapor
c Psih: Correction for heat
c z0h : Roughness length for turbulent transfer of heat [m]
c z0q : Roughness length for turbulent transfer of water vapor [m]
c
      double precision Ht,Hw,Hr,Psim,Psih,Psiq,PsimT
      double precision L,Tbar,sqrtCd,sqrtCdn
      double precision zeta,wspT
      double precision deltaB,deltaT,deltaE,D,zscale,freelim,nu
      double precision Dv,zscalev,freelimv,z0h,z0q,zetal,zetah
      double precision tstar,qstar,dumM,dumT,zetaq,zetat,y,yprime
      double precision DpsitDzetat,DpsimDzetat,con,dum
      double precision wstar,Rb,Rbw,c,qlat0,qsen0

      integer i,it
c
c Functions
c
      double precision fqlat,fqsen,dutch
      double precision unstable
c
c 1.Checks to be sure instrument heights are above snow (or grass)
c
       do 683 i=1,3
       if (height(i).le.snowdepth)then
         write(*,*) '  **Met Instruments Are Not Above Snow!!**'
         stop '**Error in Layer.in, Line 3.  Execution Terminated**'
       endif
683    continue

c
c  2. Next computes drag coefficient over snow at neutral stability
c     dlogW,dlogT and dlogQ over soil are assumed constant for the run and  
c     are computed in CALCONSTANT.f
      if(ltype .eq. 1)then
         Ht=height(1)-snowdepth
         Hw=height(2)-snowdepth
         Hr=height(3)-snowdepth
         dlogW=alog(sngl(Hw/z0))
         sqrtCdn=.4d0/dlogW
      else
         Ht=height(1)
         Hw=height(2)
         Hr=height(3)
         sqrtCdn=.4d0/dlogW
         z0h=z0
         z0q=z0
      endif
c
c
c 3. Use values at neutral stability for first guess 
c
      Psim=0d0
      Psih=0d0
      Psiq=0d0
      deltaT=Tkair-Tsurface
      deltaE=Eair-Esurface
      Tbar=(Tkair+Tsurface)/2d0
      Rbw=9.8*Ht*deltaT/(Tbar*wsp*wsp)
      zetat=Rbw*dlogW  ! A first guess
      L=Ht/zetat
c     Gas constant for dry air (287J/kg-K) is from "Handbook of 
c     Meteorology,"Berry, et al, 1945, McGraw-Hill, p. 353.
      rhoair=1d2*bpress/(287d0*tkair)
      con=100d0/(Rw*tkair*rhoair)

C
      it=0
      c=Hw/Ht
      if( Rbw .gt.1.427d0/(c*c))then !Exceeds critical Richardson number
        qsen0=1d-12
        qlat0=1d-12
        goto 60
      endif
      if(L .ge. 0d0)then  ! Bounds on zetat for N-R scheme
       zetal=0d0
       zetah=500d0
      endif
c
c 4. Compute stability functions using iterative solution.

c     Next begins Newton-Raphson iterative loop.  Loop at least twice.
 
10    Continue
      it=it+1
c 4.a Compute momentum stability functions, Cd, and wstar
      zeta=Hw/L
      IF (L .ge. 0d0)  THEN   !Stable   
        Psim=Dutch(zeta)
      ELSE  !Unstable 
        Psim=Unstable(L,Hw,z0h,1)
      ENDIF 
      sqrtCd=0.4d0/dmax1((dlogW - Psim),0.1d0)
      wstar=sqrtCd*wsp
c 4.b Estimate scalar roughness lengths and stability functions
      if(ltype .eq. 1)
     &call andreas(wstar,z0,Tkair,z0h,z0q,dlogT,dlogQ,Ht,Hr)

      zetat=Ht/L
      zetaq=Hr/L
      IF (L .ge. 0d0)  THEN   !Stable   
        if(dabs(Ht-Hw) .lt. 0.01d0)then	
           Psih=Psim
        else
           Psih=Dutch(zetat) 
        endif
        PsimT=Psih
        if(dabs(Ht-Hr) .lt. 0.01d0)then
           Psiq=Psih
        else
           Psiq=Dutch(zetaq)
        endif
      ELSE  !Unstable 
        PsimT=Unstable(L,Ht,z0h,1)
        psih=Unstable(L,Ht,z0h,2)
        if(dabs(Ht-Hr) .lt. 0.01d0)then
         psiq=psih
        else
         psiq=Unstable(L,Hr,z0h,2)
        endif
      ENDIF 
c     endif
      tstar=0.4d0*deltaT/(dlogT-Psih)
      qstar=0.4d0*deltaE*con/(dlogQ-Psiq)

c 4.c Compute Monin-Obukhov length (L)

      if(tstar .ne. 0d0 .and. wstar .ne. 0d0)then
       L=1d0 + 0.61d0*Tbar*qstar/tstar     
       dum=L
       L=L*9.8d0*0.4d0*tstar/(Tbar*wstar*wstar)
       L=1d0/L
      elseif(tstar .eq. 0d0 .and. wstar .ne. 0d0)then
       L=1d14
      endif

c 4.d Limits for extremely stable and unstable stratification 

      freelim=0d0
      freelimv=0d0
      if(zetat .gt. 500d0 .or. Ht/L .gt. 500d0)then  !very stable
        qsen0=1d-12
        qlat0=1d-12
        goto 60
      elseif(L .lt. 0 .and. L .gt. -1.0d0 )then  !very unstable
c    Compute free convection limits.  Applies generally for L < -0.1
c     Andreas and Cash, Convective heat transfer over wintertime
c     leads and polynyas, JGR 104, C11, 25,721-25,734, 1999.
c     Method is for a smooth surface with unlimited fetch.
       deltaB=(-9.8d0*deltaT/Tbar)*dum
       nu=1.315d-5 ! kinematic viscosity if air
       D=0.024/(cdryair*rhoair)
       Dv=1.1494*D
       zscale=(nu*D/deltaB)**(1d0/3d0)
       zscalev=(nu*Dv/deltaB)**(1d0/3d0)
       freelim=-0.15*rhoair*cdryair*D*deltaT/zscale
       freelimv=-0.15*Dv*deltaE*con*2.838d6/zscalev
       if(ltype .gt. 1)freelimv=(dlogT/dlogQ)*freelimv
      endif
c 4.e Compute change in zeta
      y=(Ht/L) - zetat
      if(L .gt. 0d0 .and. y .gt. 0)zetal=zetat !Update lower bound
      if(L .gt. 0d0 .and. y .le. 0)zetah=zetat !Update upper bound
      if(dabs(y) .lt. 1d-2 .and. it .gt. 2) goto 50
      if(it .gt. 50)then
        write(80,*) 'Convergence problem in QTURB.f'      
        goto 50
      endif
c 4.f Use Newton-Raphson scheme for stable,sequential estimates for
c     unstable
      if(L .ge. 0d0)then
c       Next is derivitive of Psi
        DpsitDzetat=
     &     -0.70-0.75d0*dexp(-0.35d0*zetat)*(6d0-0.35d0*zetat)
        DpsimDzetat= 
     &     c*(-0.70-0.75d0*dexp(-0.35d0*zeta)*(6d0-0.35d0*zeta))
        dumM=1d0/(dlogW-Psim)
        dumT=1d0/(dlogt-Psih)
c       Next is derivitive of y (ignoring humidity effects)
        yprime=(Ht/L)*(DpsitDzetat*dumT-2d0*DpsimDzetat*dumM) - 1d0
        zetat=zetat - (y/yprime)
        if(zetat .le. zetal .or. zetat .ge. zetah)then
c       Root out of range. Bisection instead of Newton-Raphson scheme
         zetat=(zetah+zetal)/2d0
        endif
      else
        zetat = Ht/L
      endif
      L=Ht/zetat
      goto 10

c
c 5.  Final computation of qsen and qsen. Fluxes are computed in MAIN.
c Note qlat and qsen add the windless coefficients and qlat is times
c the latent heat of sublimation or evaporation.
c
50    continue
      qlat0=fqlat(dlogW,dlogQ,Rw,tkair,Psim,Psiq,wsp)
      qsen0=fqsen(dlogW,dlogT,cdryair,rhoair,Psim,Psih,wsp)
60    if(dliqvol .gt. 0.020 .or. tsurface .gt. 273.15)then
c       For wet snow, use latent heat of evaporation
        qlat= Ck +qlat0*2.5045d6
      else
c       For dry snow, use latent heat of sublimation
        qlat=Ck +qlat0*2.838d6
      endif
      if(L .ge. 0d0)then
        qsen=qsen0  ! March 19,2002
        if(qsen0 .lt. dmax1(Csk,1d-1) .and. qsen0 .gt. 0)
     &    qsen = dmax1(Csk,1d-1)
      else
        dum=0d0
        if(deltaT .ne. 0d0)dum=-freelim/deltaT
        qsen=dmax1(qsen0,dum)
        dum=0d0
        if(deltaE .ne. 0d0)dum=-freelimv/deltaE
        qlat=dmax1(qlat,dum)
      if(L .lt. 0 .and. L .gt. -1.0d0 )then  !very unstable
       endif
      endif
c
c 6. Compute bulk Richardson number
c
c     First adjust windspeed to height of temperature measurement
      if( L .gt. 0d0)then
       PsimT=Psih
      else
       PsimT=Unstable(L,Ht,z0h,1)
      endif
      wspT=(wsp/(dlogW-Psim))*(dlog(Ht/z0)-PsimT)
      Rb=9.8*Ht*deltaT/(Tbar*wspT*wspT)

      return
      end

c***********************************************************************
c Next are function subroutines
c***********************************************************************

c Function FQSEN
      double precision function fqsen(dlogW,dlogT,cdryair,rhoair,Psim,
     &  Psih,wsp)
      double precision dlogW,dlogT,cdryair,rhoair,psim,psih,wsp
      fqsen = cdryair*rhoair*wsp*0.16d0/((dlogW - Psim)*(dlogT - Psih))
      return
      end  

c Function FQLAT
      double precision function fqlat(dlogW,dlogQ,Rw,Tkair,psim,psiq
     &  ,wsp)
      double precision dlogW,dlogQ,Rw,Tkair,psim,psiq,wsp
c     Since vapor pressure is in mb, fqlat is multiplied by 100 to
c     convert to Pascals
      fqlat = wsp*0.16d2/(Rw*tkair*(dlogW - Psim)*(dlogQ - Psiq))
      return
      end 

c Function DUTCH
      double precision function Dutch(zeta)
      double precision zeta
      Dutch=-(0.70d0*zeta+0.75d0*(zeta-14.28d0)*	
     &       dexp(-0.35d0*zeta)  + 10.71d0)	
      return
      end 

c Function UNSTABLE
      double precision function Unstable(L,Ht,zt,icall)
      double precision zeta,x,L,Ht,Llim,zt
      integer icall      
c Set an upper limit on function at L = -.1 or  at 200 times
c the large -Ht/L limit where log(-4*Ht/L) is 0.  This asymptotic
c limit is from Godfrey and Beljaars (1991).  This is a semi-
c empirical fix to keep the heat flux from increasing with 
c decreasing -L.
      Llim=dmin1(-1d-1,-100d0*zt)
      zeta=Ht/dmin1(L,Llim)
      x=(1d0-16d0*zeta)**0.25d0
      if(icall .le. 1)then      
       Unstable=dlog((1+x*x)/2d0)+2d0*dlog((1d0+x)/2d0)- 2d0*datan(x)
     &          +1.5707963d0
      else      
       Unstable=2d0*dlog((1+x*x)/2d0)
      endif
      return
      end      
c***********************************************************************
      subroutine andreas(wstar,z0,Tkair,z0h,z0q,dlogT,dlogQ,Ht,Hr)

c Arguments
      double precision wstar,z0,Tkair,z0h,z0q,dlogT,dlogQ,Ht,Hr
c Local
      double precision Reynolds,b0,b1,b2,dlogReynolds,dum,dlogsqrd

c    Compute scaler roughness lengths using procedure
c    in Andreas, 1987, Boundary-Layer Meteorology 38, 159-184.  Only
c    use this procedure for snow.

c     Expression for kinematic viscosity of air is taken from program
c     of Launiainen and Vihma, 1990, Environmental Software, vol. 5,
c     No. 3, pp. 113 - 124.

c Heat
      Reynolds=z0*wstar*1d7/(.9065*Tkair-112.7d0)
      if(Reynolds .le. 0.135d0)then
       b0=1.250
       b1=0d0
       b2=0d0
      elseif(Reynolds .lt. 2.5d0)then
       b0=0.149d0
       b1=-0.550d0
       b2=0d0
      else
       b0=0.317d0
       b1=-0.565d0
       b2=-0.183d0
      endif
      dlogReynolds=dlog(Reynolds)
      dlogsqrd=dlogReynolds*dlogReynolds
      dum=b0 + b1*dlogReynolds + b2*dlogsqrd
      z0h=z0*dexp(dum)

c Moisture 
      if(Reynolds .le. 0.135d0)then
       b0=1.61d0
       b1=0d0
       b2=0d0
      elseif(Reynolds .lt. 2.5d0)then
       b0=0.351d0
       b1=-0.628d0
       b2=0d0
      else
       b0=0.396d0
       b1=-0.512d0
       b2=-0.180d0
      endif
      dum=b0 + b1*dlogReynolds + b2*dlogsqrd
      z0q=z0*dexp(dum)

      dlogT= dlog(Ht/z0h)
      dlogQ= dlog(Hr/z0q)

      return
      end
