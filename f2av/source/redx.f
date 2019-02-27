      REAL FUNCTION REDX(X,Q2,ITYPE,CME,ICUD)
C------------------------------------------------------------------------
C IYTPE:    1  NC e-
C           2  CC e-
C           3  NC e+
C           4  CC e+
C
C 
C CME:     820 or 920
C          
C
C ICUD:     1  up
C           0  central
C          -1  down 
C
C
C SG: Add Fractal fit for low Q2
C
C------------------------------------------------------------------------
      implicit none
      INTEGER I,J,Ifile,ITYPE,ICUD,icme
      real cme
      REAL    ZPDFN
C                         IREA, CME 
      REAL    PDFG(161,161,4  ,2   )
      REAL    PDFU(161,161,4  ,2   )
      REAL    PDFD(161,161,4  ,2   )
      REAL    X,Q2
 
      COMMON/CMNPDG/PDFG 
      LOGICAL LFIRST
      DATA    LFIRST /.TRUE./
      SAVE    LFIRST
      
      character*12 CFILES(4,2)
      data CFILES/
     $     'emnc_820.dat','emcc_820.dat','epnc_820.dat','epcc_820.dat',
     $     'emnc_920.dat','emcc_920.dat','epnc_920.dat','epcc_920.dat'/
      
      character*44 CDIR1
      parameter (CDIR1='/afs/desy.de/group/h1zeus/www/html/inclusive')
      character*26 CDIR2
      parameter (CDIR2='/internal/code/herapdf0.2/')
      character*132 FNAM
C------------------------------------------------------------------------
      double precision pi, alem,conv, xmw, gmuer
      parameter (pi=3.14159267,alem=1.0/137.0359895)
      parameter (conv=0.3893834e+09,xmw=80.41,gmuer=1.16807E-05)
      
C Below Q2fract1: pure fractal, above Q2Fract2: pure HERAPDF0.2,
C                 linear param. in between
      real Q2Fract1,Q2Fract2
      
      parameter (Q2Fract1=2.5,Q2Fract2=3.0)
      real y,fact,fract

      real red920, red820, y920, y820, yp920, yp820, yp
      real fitfractH
C------------------------------------------------------------------------

C     Read in PDF values
      IF (LFIRST) THEN
         do icme=1,2
            do ifile=1,4
             FNAM = CDIR1//CDIR2//CFILES(ifile,icme)
             OPEN(UNIT=50,FILE=FNAM,STATUS='OLD')
             print '(''Read from file='',A90)',FNAM
             DO I=1,161
                DO J=1,161
                   READ(50,*) PDFG(J,I,ifile,icme)
     $                  ,PDFU(J,I,ifile,icme),PDFD(J,I,ifile,icme)
                ENDDO
             ENDDO
             close (50)
C             print *,PDFG(5,1,ifile,icme),ifile,icme
          enddo
        enddo
        LFIRST = .FALSE.
      ENDIF
      

      if (CME.gt.260.) then

         if (CME.gt.310.) then
            ICME = 2
         else
            ICME = 1
         endif

C      print *,'GOTCME',CME,ICME

         IF (ICUD.EQ.0) THEN
            REDX = ZPDFN(X,Q2,PDFG(1,1,Itype,ICME))
         ELSEIF (ICUD.EQ.1) THEN
            REDX = ZPDFN(X,Q2,PDFU(1,1,Itype,ICME))
         ELSEIF (ICUD.EQ.-1) THEN
            REDX = ZPDFN(X,Q2,PDFD(1,1,Itype,ICME))
         ELSE
            PRINT*, 'WRONG ICUD Option'
            STOP
         ENDIF
      else
C
C Low energy run
C
         red920 = ZPDFN(X,Q2,PDFG(1,1,Itype,2)) 
         red820 = ZPDFN(X,Q2,PDFG(1,1,Itype,1)) 
         y920 = q2/(4.*920.*27.5*x)
         y820 = q2/(4.*820.*27.5*x)
         y    = q2/(CME**2*x)
         yp920 = y920**2/(1+(1-y920)**2)
         yp820 = y820**2/(1+(1-y820)**2)
         yp    = y**2/(1+(1-y)**2)

         REDX = red920 + (red920-red820)/(yp820-yp920)*(yp920-yp)

C         print *,'hoho',red920,red820,redx,yp,yp920,yp820
      endif
      
      if (IType.eq.2 .or. IType.eq.4) then      
         fact = xmw**4/(xmw**2+q2)**2*GMUER**2/(2*pi*x)
         redx = redx * fact * conv
      endif

      if (IType.eq.3) then
         if (Q2.lt.Q2fract2) then
            Y = q2/(CME**2*X)
            fract = fitfractH(x,q2,y)
            if (Q2.le.Q2Fract1) then
               redx = fract
            else
               redx = (redx*(q2-q2fract1)+fract*(q2fract2-q2) )
     $              /(q2fract2-q2fract1)
            endif
         endif
      endif

C      RETURN
      END
C------------------------------------------------------------------------

      real Function fitfractH(xv,q2v,yv)
      implicit none
      real xv,q2v,yv
      logical lfirst 
      data lfirst/.true./
      double precision pval(10)
      double precision D0,D1,D2,D3,Q0,R
      double precision F2,FL,XS
      integer i
      character*44 CDIR1
      parameter (CDIR1='/afs/desy.de/group/h1zeus/www/html/inclusive')
      character*26 CDIR2
      parameter (CDIR2='/internal/code/herapdf0.2/')
      character*132 FNAM
C---------------------------------------
      if (lfirst) then
         lfirst = .false.
         print *,'Init fit fract'
         FNAM = CDIR1//CDIR2//'fractout.dat'
         open (33,file=FNAM,status='old')
         read (33,'(6E16.8)') (pval(i),i=1,6)
         close (33)
         D0 = pval(1)
         D1 = pval(2)
         D2 = pval(3)
         D3 = pval(4)
         Q0 = pval(5)
         R  = pval(6)
         print *,'Fractal parameters:'
         print 18,(pval(i),i=1,6)
 18      format(6F12.6)
      endif

      if (yv.gt.1.) then
C         fitfract = 0.
C         Return
      endif

      F2=D0*Q0*((Q2v/(Q2v+Q0))**(D2-1))
     $     *(xv**(-D2+1))/(1+D3-D1*log(xv))
      F2=F2*(xv**(-D1*log(1+Q2v/Q0))*((1+Q2v/Q0)**(D3+1))-1)
      
      FL = F2*R/(1+R)
      
      XS = F2 - Yv**2/(1+(1-Yv)**2) * FL
      
C      print *,xs,pval

      FitFractH = XS

      end
