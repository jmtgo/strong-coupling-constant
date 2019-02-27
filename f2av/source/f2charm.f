*
*      GET F2C from RIEMERSMA
*
       FUNCTION F2CHAR(Q2CEN,XCEN,SCMS)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       REAL Q2CEN,XCEN,SCMS
       parameter (pi=3.14159267D0,alem=1.0D0/137.0359895D0)
       parameter (conv=0.3893834e+09)

* select the PDF 1=MRST FF 2=CTEQ5F3
       COMMON / USEPDF/ IPDF
       INTEGER IPDF
       INTEGER IQUARK
       COMMON / GROUP / CA, CF, TF 
* steering 1/2 c/b
      IQUARK=1
* select the PDF 1=MRST FF 2=CTEQ5F3
      IPDF=1
*      IPDF=2

* CTEQ requires an initialisation
      IF (IPDF.EQ.2) THEN
         CALL SETCTQ5(6)
      ENDIF

*
* ..QCD colour factors
*
      CA = 3.D0
      CF = 4.D0/3.D0
      TF = 1.D0/2.D0
*
* ..Charge, mass, Q^2 and factorization scale (all fixed here)
*
      IF (IQUARK.EQ.1) THEN
         EC = 2.D0/3.D0
* cteq use MC=1.3 MRST use mc=1.43
         IF (IPDF.EQ.1) THEN
            MC = 1.43D0
         ELSE
            MC = 1.3D0
         ENDIF
         
      ELSE
         EC = 1.D0/3.D0
         IF (IPDF.EQ.1) THEN
* mrst use mb=4.3
            MC = 4.3D0
         ELSE
* cteq use mb=4.5
            MC = 4.5D0
         ENDIF 
      ENDIF
      
*
      AM = 4.D0
      AQ = 1.D0

*     
C      write(*,*)'F2CHARM ist gerufen!',Q2CEN,XCEN,SCMS  
      q2c=Q2CEN
      xc=XCEN
      
      y=q2c/scms/xc
      yplus =1.0+(1.0-y)**2
      yminus=1.0-(1.0-y)**2
      fac = xc*q2c*q2c/(2.0*pi*alem*alem*yplus*conv)
      MF = SQRT (Q2C+4*MC*MC)
      MR = MF

*     Calculate F2c 
C      write(6,*),'fctn:scales',MF,MR,MC,EC
C      write(6,*),'Q2,x,cms',q2c,xc,scms
      F2CHAR = F2H(XC, Q2C, MF, MR, MC, EC)
C      write(6,*),'fctn:',q2c,xc,scms,f2char
      RETURN
      END
