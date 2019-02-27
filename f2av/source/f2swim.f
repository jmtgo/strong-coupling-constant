      Subroutine F2SWIM(Reaction,Qin,Xin,Qout,Xout,Scale,YMAXe
     $     ,POLAR,CME)
C-----------------------------------------------------------------------
C
C Created 5 July 2004 by SG.
C Determine the scale factor needed to rescale 
C
C   Xsect (Qin,Xin) --> Xsect (Qout,Xout)
C
C Uses XS parametrization as obtained in H1 QCD fit ! 
C
C Modification:
C   25 June 06 -- add correction of single diff. cross section for upper Y
C              limit. (YMaxe)
C   28 June 06 -- add polarization
C 
C-----------------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      integer Reaction
      real Qin,Xin,Qout,Xout,Scale,YMaxE,Polar,CME
      real F2a,F2b,ToPolar,Yin,Yout
C
C -- H1 parametrizations of NC/CC cross sections:
C
      real signc,sigcc,signc2,f2fract
C
C --- From Eram parameterization of integrated dSdQ2
C
      real dsdq2,FitFract
      integer irea,IRsav

      real F2char
C HERA0.2 parameterization:
      real redx
C---------------------------------------------------------------

C Save IRSTUDY flag:
      IRSav = IRSTUDY
C Switch off IRSTUDY:
      IRSTUDY = 0

C 5 July 2004, SG: for now no swiming: 
      scale = 1.0
 
      Yin  = Qin /(CME*Xin)
      Yout = Qout/(CME*Xout)

C 28 June 2006: Find desired polarization for the measurement:
      ToPolar = 0.
      do irea=1,NREACT
         if (Reaction.eq.RLIST(irea)) then
            ToPolar = RPOLAR(irea)
         endif
      enddo

      if (abs(Qin-Qout).gt.0.001.or.abs(Xin-Xout).gt.1.E-7
     $     .or.YMAXINT.ne.YMAXe
     $     .or.(ToPolar.ne.Polar)) then

         if (Reaction.eq.415) then
C NC CHARM PRODUCTION: SWIMMING via RIEMERSMA, CTEQ5FF3, Mc=1.3
C            write(6,*)'swimming: Q2, x CMS start',qin,xin,CME
C            write(6,*)'swimming: Q2, x CMS end',qout,xout,CME
            F2a = F2CHAR(qin,xin,CME)
            F2b = F2CHAR(qout,xout,CME)
C            write(6,*)'swimming via riemersma',F2a,F2b

         elseif (Reaction.eq.515) then
C NC e- p

            
            if (LHERA02) then
               F2a = redx(xin,qin,1,sqrt(CME),0)
               F2b = redx(xout,qout,1,sqrt(CME),0)
            else
               F2a = signc2(xin,qin,-1.,CME)
               F2b = signc2(xout,qout,-1.,CME)
            endif
         elseif (Reaction.eq.615) then
C NC e+ p


            if (LFractal) then
               F2a = FitFract(xin,qin,yin)
               F2b = FitFract(xout,qout,yout)
            else
C ADD LOW Q2 FRACTAL FIT
               if (LHERA02) then
C Fractal is already in:
                  F2a = redx(xin,qin,3,sqrt(CME),0)
                  F2b = redx(xout,qout,3,sqrt(CME),0)                     
               else
                  if (qin.ge.4) then
                     F2a = signc2(xin,qin,+1.,CME)
                     F2b = signc2(xout,qout,+1.,CME)
                  else
                     F2a = FitFract(xin,qin,yin)
                     F2b = FitFract(xout,qout,yout)
                  endif
               endif
            endif

         elseif (Reaction.eq.3515) then
C CC e- p
            if (LHERA02) then
               F2a = redx(xin,qin,2,sqrt(CME),0)
               F2b = redx(xout,qout,2,sqrt(CME),0)                
            else
               F2a = sigcc(xin,qin,-1.,CME)
               F2b = sigcc(xout,qout,-1.,CME)
            endif
         elseif (Reaction.eq.3615) then
C CC e+ p
            if (LHERA02) then
               F2a = redx(xin,qin,4,sqrt(CME),0)
               F2b = redx(xout,qout,4,sqrt(CME),0)
            else
               F2a = sigcc(xin,qin,+1.,CME)
               F2b = sigcc(xout,qout,+1.,CME)
            endif

C 22 June 06: Add polarized integrated cross sections
         elseif (Reaction.eq.10615) then
C NC e+ L polarization:            
            if (xin.eq.3.0) then ! Flag that Q2 only
               F2a = dsdq2(2,qin,YMaxE,0.0,0.,920.,27.6,+1.,POLAR)
               F2b = dsdq2(2,qout,YMaxInt,0.0,0.,920.,27.6,+1.,ToPOLAR)
            else
               F2a = 1.
               F2b = 1.
            endif
         elseif (Reaction.eq.20615) then
C NC e+ R polarization:            
            if (xin.eq.3.0) then ! Flag that Q2 only
               F2a = dsdq2(2,qin,YMaxE,0.0,0.,920.,27.6,+1.,POLAR)
               F2b = dsdq2(2,qout,YMaxInt,0.0,0.,920.,27.6,+1.,ToPOLAR)
            else
               F2a = 1.
               F2b = 1.
            endif
         elseif (Reaction.eq.10515) then
C NC e- L polarization:            
            if (xin.eq.3.0) then ! Flag that Q2 only
               F2a = dsdq2(2,qin,YMaxE,0.0,0.,920.,27.6,-1.,POLAR)
               F2b = dsdq2(2,qout,YMaxInt,0.0,0.,920.,27.6,-1.,ToPOLAR)
            else
               F2a = 1.
               F2b = 1.
            endif
         elseif (Reaction.eq.20515) then
C NC e+ R polarization:            
            if (xin.eq.3.0) then ! Flag that Q2 only
               F2a = dsdq2(2,qin,YMaxE,0.0,0.,920.,27.6,-1.,POLAR)
               F2b = dsdq2(2,qout,YMaxInt,0.0,0.,920.,27.6,-1.,ToPOLAR)
            else
               F2a = 1.
               F2b = 1.
            endif
         else
            print *,'UNKNOWN REACTION=',REACTION
            print *,'NO INTERPOLATION'
            F2a = 1.
            F2b = 1.
         endif



         Scale = F2b/F2a

         print *,'Interpolation for reaction =',REACTION
         print '(2F13.7,'' -> '',2F13.7,'' Scale='',F10.4)',
     $        Qin,Xin,Qout,Xout,Scale

         if (YmaxInt.ne.YmaxE) then
            print
     $ '(''Correct for Y integration range'',F5.2,'' ->'',F5.2)',
     $           YmaxE,YmaxInt
         endif
         if (Polar.ne.ToPolar) then
            print
     $ '(''Correct for Polarization'',F8.4,'' ->'',F8.4)',
     $           Polar,ToPolar
         endif

      endif
C Restore IRSTUDY:
      IRSTUDY = IRSav
      end



Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      real function signc(x,q2,qin,S)
*.....***********************************
      parameter (pi=3.14159267,alem=1.0/137.0359895)
      parameter (conv=0.3893834e+09,xmw=80.41,gmuer=1.16807E-05)
      real scms,qin
*
*      qin (real)   0.0/+-1.0   charge of lepton . Taken from
*                           vector qlep if set to 0.0
*
*

      ypol = 0.  ! no po;arization

 
*.....for H1 2000 PDF fit polarised or not!
      redxsec = signcpol(x,q2,qin,ypol,S)

      signc = redxsec

      end

C---------------------------------------------
      real function signc2(x,q2,qin,s)
C
      implicit none
      include '../include/f2.inc'
      
      real x,q2,qin,s
      real y,Yp,Ym
      real f2,fl,xf3
C H1 parametrizations:
      real f2nchiq200,flnchiq200,f3nchiq200
      real signc
C-------------------------------------------------------------
      if (q2.le.4.) then
         signc2 = signc(x,q2,qin,s)
         Return
      endif
      y = q2/s/x
      y = min(y,0.99)
      yp = 1 + (1-y)**2
      ym = 1 - (1-y)**2

C
      f2 = f2nchiq200(x,q2)
      fl = flnchiq200(x,q2)
      xf3 = f3nchiq200(x,q2)

      if (IRSTUDY.eq.1) then
         FL = 0
      endif

CC      print *,x,q2,y,f2,fl,xf3
C
      if (qin.gt.0) then
         signc2 = f2 - Ym/Yp * xF3 - y*y/Yp * FL 
      else
         signc2 = f2 + Ym/Yp * xF3 - y*y/Yp * FL 
      endif
C-------------------------------------------------------------
      end

      real function sigcc(x,q2,qin,S)
C
C 
C
      parameter (pi=3.14159267,alem=1.0/137.0359895)
      parameter (conv=0.3893834e+09,xmw=80.41,gmuer=1.16807E-05)
C-----------------------
      sigcc=0.0

      scms = S
 
      y      = q2/scms/x
      yplus  = 1.0+(1.0-y)**2
      yminus = 1.0-(1.0-y)**2


 
      if (qin.gt.0.0) then
         w2 = w2cphiq200(x,q2)
         wl = wlcphiq200(x,q2)
         xw3= w3cphiq200(x,q2)   
      else
         w2 = w2cehiq200(x,q2)
         wl = wlcehiq200(x,q2)
         xw3= w3cehiq200(x,q2)         
      endif   
      
* ER  Note there is an extra factor of 2 in the definition of the
*     reduced cross-section, this is because the structure functions are
*     defined using the hector prescription where e.g. W2=2d0*(db+sb+u+c)
*     (routine STRUFC). That definition is used to create the
*     parametrisation functions of the S.F. in QCDNUM, and so MUST be
*     adopted here.

      pol = 0.
      sfterm  = (yplus*W2  - qin*yminus*XW3  - y**2*WL )
 
      dsdxdq2  = sfterm   * conv * (1.0+qin*pol)/2.0
      redxs = dsdxdq2 / (2.0  * conv )
 
      sigcc=redxs

C 21 June: Data is NOT reduced CC:

      fact = xmw**4/(xmw**2+q2)**2*GMUER**2/(2*pi*x)

      sigcc = sigcc * fact * conv
 
      return
      end

      real function signcpol(x,q2,ech,polar,S)
      implicit none
      real x,q2,ech,polar,y,yp,ym,s
      real sinehp00h1ms,sinphp00h1ms
      real sinehn00h1ms,sinphn00h1ms
      logical first
      data first/.true./
 
      if(first) then
        first=.false.
      endif
 
      y=q2/(x*s)
      if(y.lt.0.or.y.gt.1.) then
         signcpol=0.
         return
      else
         yp=1.+(1.-y)**2
         ym=1.-(1.-y)**2
      endif
 
      if(ech.lt.0.) then
        signcpol=(1.-polar)/2.*sinehn00h1ms(x,q2)+
     &           (1.+polar)/2.*sinehp00h1ms(x,q2)
      else
        signcpol=(1.-polar)/2.*sinphn00h1ms(x,q2)+
     &           (1.+polar)/2.*sinphp00h1ms(x,q2)
      endif
      return
      end

       function f2fract(X,Q2,FLAG)
       real X,xx,Q2,F2
       real D0,D1,D2,D3,Q0
       integer FLAG
       if (flag.le.2) then
          if (flag.eq.1) then
* all param fit
             D0 = 1.40403e+00
             D1 = 7.30814e-02
             D2 = 1.01308e+00
             D3 =-1.28739e+00
             Q0 = 6.19745e-02
          else
* D2 fixed to 1
             D0 = 1.68745e+00
             D1 = 7.39626e-02
             D2 = 1.00000e+00
             D3 =-1.28158e+00
             Q0 = 5.10077e-02
          endif
          F2=D0*Q0*(x**(-D2+1))/(1+D3-D1*log(x))
          F2=F2*(x**(-D1*log(1+Q2/Q0))*((1+Q2/Q0)**(D3+1))-1)
       else if (flag.le.4) then
          if (flag.eq.3) then
*
             D0 = 6.42081e-01
             D1 = 6.35265e-02
             D2 = 0.075
             D3 =-1.28525e+00
             Q0 = 1.35189e-01
             F2=D0*Q0*((Q2/(Q2+Q0))**D2)*(x**(-D2))/(1+D3-D1*log(x))
             F2=F2*(x**(-D1*log(1+Q2/Q0))*((1+Q2/Q0)**(D3+1))-1)
          else 
*
             D0 = 1.06700e+00
             D1 = 6.70924e-02
             D2 = 7.69720e+01
             D3 =-1.40026e+00
             Q0 = 1.16272e-01
             F2=D0*Q0*((Q2/(Q2+D2))**0.07)*(x**(-0.07))/(1+D3-D1*log(x))
             F2=F2*(x**(-D1*log(1+Q2/Q0))*((1+Q2/Q0)**(D3+1))-1)
          endif
       else
          if (flag.eq.5) then
* x modified, fit to BPT + H1, no H1 renormalisation
*             D0 = 5.99407e-01
*             D1 = 6.41714e-02
*             D2 = 1.07967e+00
*             D3 =-1.29723e+00
*             Q0 = 1.47090e-01
*             rm = 1.78065e-01
* another fit
*             D0 = 5.70948e-01
*             D1 = 7.13773e-02
*             D2 = 1.07008e+00
*             D3 =-1.36720e+00
*             Q0 = 1.90933e-01
*             rm = 1.67041e-01
*
* last fit, gives about the same as the first one
*
             D0 = 5.99405e-01
             D1 = 6.41721e-02
             D2 = 1.07967e+00 
             D3 =-1.29723e+00
             Q0 = 1.47093e-01
             RM = 3.56117e-01
*
             RMP2 = 0.938*0.938
             RW2  = RMP2+Q2*(1-x)/x
c             xx = x+4*rm*rm/(RW2-RMP2)
             xx = x + rm*rm/(RW2-RMP2)
             F2=D0*Q0*(xx**(-D2+1))/(1+D3-D1*log(xx))
             F2=F2*(xx**(-D1*log(1+Q2/Q0))*((1+Q2/Q0)**(D3+1))-1)
          else
             F2=1.
          endif
       endif
       f2fract=F2
       end

      real Function fitfract(xv,q2v,yv)
      implicit none
      include '../include/f2.inc'
      real xv,q2v,yv
      logical lfirst 
      data lfirst/.true./
      double precision pval(10)
      double precision D0,D1,D2,D3,Q0,R
      double precision F2,FL,XS
      integer i
C---------------------------------------
      if (lfirst) then
         lfirst = .false.
         print *,'Init fit fract'
         D0 = 0.67558497E+00
         D1 = 0.55442097E-01  
         D2 = 0.10800000E+01
         D3 =-0.11988581E+01  
         Q0 = 0.11159892E+00 
         R  = 0.58655491E+00 
      endif

      if (yv.gt.1.) then
C         fitfract = 0.
C         Return
      endif

      F2=D0*Q0*((Q2v/(Q2v+Q0))**(D2-1))
     $     *(xv**(-D2+1))/(1+D3-D1*log(xv))
      F2=F2*(xv**(-D1*log(1+Q2v/Q0))*((1+Q2v/Q0)**(D3+1))-1)
      

C
C Add code to study influence of R assumption on average:
C
      if (IRSTUDY.eq.0) then
         FL = F2*R/(1+R)
      else if (IRSTUDY.eq.1) then
         FL = 0
      else if (IRSTUDY.eq.2) then
         FL = 0.5*F2     ! assume R=1
      else
         print *,'FitFract: wrong IRSTUDY --- stop'
         stop
      endif
      
      XS = F2 - Yv**2/(1+(1-Yv)**2) * FL
      
C      print *,xs,pval

      FitFract = XS

      end
