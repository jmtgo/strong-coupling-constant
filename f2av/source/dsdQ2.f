*.....****************************************************************
      real function dsdQ2(icc,q2c,ymax,Emin,Ptmin,Eprot,Elep,Qlep,Pol)
*.....****************************************************************

      implicit none
      
*     integration for dsigma/dq2 NC or CC

*     icc: 1 = CC
*          2 = NC
*     ymax  cut applied for NC and CC
*     Emin  cut applied for NC only
*     Ptmin cut applied for CC only
*     Eprot / Elep = proton/lepton beam energies
*     Qlep = lepton charge +/-1
*     Pol  = lepton polarisation -1 < Pol < +1

*     theta max cut implemented but not used for NC


      real scms,sigtot,q2min,q2max,xmin,xmax,ymax,emin,thmin,ptmin
      real dq2,dx,q2width,xwidth,Qlep,ptc
      real Pol,Eprot,Elep

      real xc,q2c,ec,thc,yc
      real yplus,yminus,fac
      real pi,alem,conv,xmw,gmuer
      parameter (pi=3.14159267,alem=1.0/137.0359895)
      parameter (conv=0.3893834e+09,xmw=80.41,gmuer=1.16807E-05)

      real xsec,signcpol,sigcc,sigtit,sigf2
      real bxmin,bxmax,ybmax
      integer nn,i,j,icc

*     scms = centre of mass energy
      scms   = 4.0*Elep*Eprot

*.....500 bins is a good number, accurate to 1.10^-4
      nn=500
      xmin=q2c/scms
      xmax=0.95

      dx =(log(xmax )-log(xmin ))/float(nn)
      dsdq2=0.0
      bxmin=1.0
      bxmax=-1.0
      do 20 j=1,nn
         xc  = log(xmin )+(float(j)-0.5)*dx   
         xc  = exp(xc)
         yc=q2c/(scms*xc)
         if (yc.gt.ymax)       goto 20
         if (icc.eq.1) then
            ptc  = sqrt( q2c*(1.0-yc) )
            if (ptc.lt.ptmin)  goto 20
            xsec = sigcc(xc,q2c,Qlep,scms)
         else
            ec=Elep*(1.0-yc)+Eprot*xc*yc
            thc=2.0*acos(sqrt(q2c/(4.0*Elep*ec)))
c           if (thc*57.3.gt.thmax) goto 20
            if (ec.lt.emin)       goto 20

c...........Calculate conversion factors from red-x-sec to ds/dxdQ2
            yplus  = 1.0+(1.0-yc)**2
            yminus = 1.0-(1.0-yc)**2
            fac = xc*q2c*q2c/(2.0*pi*alem*alem*yplus*conv)
            
            xsec = signcpol(xc,q2c,Qlep,Pol,scms)/fac

         endif
         dsdq2=dsdq2+xsec*xc*dx
         bxmin=min(bxmin,xc)
         bxmax=max(bxmax,xc)
 20   continue

      return
      end


