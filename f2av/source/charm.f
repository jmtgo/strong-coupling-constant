*
* ..Main routine for checking the LO and NLO heavy-quark F_2 and F_L
*
c       PROGRAM CHARM
       SUBROUTINE  CHARM
*
       IMPLICIT DOUBLE PRECISION (A - Z)
*
       PARAMETER ( NXX = 20 )
       DIMENSION XX(NXX)
*
       COMMON / GROUP / CA, CF, TF 
*
* ..Bjorken-x list
*
       DATA XX / 1.D-7,  5.D-7,  1.D-6,  5.D-6,  1.D-5,  5.D-5,
     1           1.D-4,  5.D-4,  1.D-3,  5.D-3,  1.D-2,  5.D-2,
     2           1.D-1, 1.5D-1,  2.D-1,  3.D-1,  4.D-1, 5.D-1,
     3           7.D-1,  9.D-1 /
*
* ..QCD colour factors
*
       CA = 3.D0
       CF = 4.D0/3.D0
       TF = 1.D0/2.D0
*
* ..Charge, mass, Q^2 and factorization scale (all fixed here)
*
       EC = 2.D0/3.D0
C       MC = 1.414213562D0
       MC = 1.5D0
*
       QS = 10.D00
       AM = 4.D0
       AQ = 1.D0
       MF = SQRT (AM* MC**2 + AQ* QS)

*
* ..Bjorken-x loop
*

       DO IXX = 1,NXX
*
          XB  = XX (IXX)

* ..Calculate F2c and FLc
          MR=MF
          F2C = F2H(XB, QS, MF, MR, MC, EC)
          FLC = FLH(XB, QS, MF, MR, MC, EC)
*
* ..Output
*
          WRITE (6,10) XB, QS, F2C, FLC
 10       FORMAT (3X,1(1PE14.5),2X,1(1PE14.5),2X,1(1PE14.5),2X
     &         ,1(1PE14.5))
*
       ENDDO
*
       STOP
       END
*
* =====================================================================
*
*
* ..Lowest-order heavy-quark contribution to F_2 (photon-gluon fusion)
*
       FUNCTION F2H (XB, QS, MF, MR, MH, EH)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
*
       PARAMETER ( PI = 3.1415926536D+00 )
       INTEGER CNT
       EXTERNAL F2HI,F2LI
*
* ..Common block passing variables to the convolution integrand
*
       COMMON / PARF2C / X, Q2, MF2, MR2, MH2, EH2
*
       X  = XB
       Q2 = QS
       MF2 = MF * MF
       MR2 = MR * MR
       MH2 = MH * MH
       EH2 = EH * EH

*
* ..The integration, if the threshold energy is exceeded 
*
       A  = 1.D0 + 4.D0 * MH2 / Q2
       AX = A * X
       IF ( AX .GE. 1.D0 ) THEN
         F2H = 0.0D0
       ELSE
         F2H = EH2 * DAIND (AX,1.D0,F2HI,1.D-5,2,10000,CNT,ERR) 
       END IF

       RETURN
       END
*
* =====================================================================
* =====================================================================
*
*
* ..Lowest-order heavy-quark contribution to F_L (photon-gluon fusion)
*
       FUNCTION FLH (XB, QS, MF, MR, MH, EH)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
*
       PARAMETER ( PI = 3.1415926536D+00 )
       INTEGER CNT
       EXTERNAL F2HI,F2LI
*
* ..Common block passing variables to the convolution integrand
*
       COMMON / PARF2C / X, Q2, MF2, MR2, MH2, EH2
*
       X  = XB
       Q2 = QS
       MF2 = MF * MF
       MR2 = MR * MR
       MH2 = MH * MH
       EH2 = EH * EH
*
* ..The integration, if the threshold energy is exceeded 
*
       A  = 1.D0 + 4.D0 * MH2 / Q2
       AX = A * X
       IF ( AX .GE. 1.D0 ) THEN
         FLH = 0.0D0
       ELSE
         FLH = EH2 * DAIND (AX,1.D0,F2LI,1.D-5,2,10000,CNT,ERR) 
       END IF

       RETURN
       END
*
* =====================================================================
*
*
* ..The convolution integrand
*
       FUNCTION F2HI (E)
*
       IMPLICIT DOUBLE PRECISION (A - H, L - Z)
*
* ..Input common-block

       PARAMETER ( PI = 3.1415926536D+00 )
*      
       COMMON / PARF2C / X, Q2, MF2, MR2, MH2, EH2

       COMMON / USEPDF/ IPDF
       INTEGER IPDF
*
* ..The parton densities at x = E (GL = x*g, as usual)

c       ISET = 2
c       CALL GRV98PA (ISET, E, MF2, MR2, UV, DV, US, DS, SS, GL)
c       write(6,*) 'ipdf',ipdf
       IF (IPDF.EQ.1) THEN
          FB  = 0D0
          FC  = 0D0
          CALL MRST2004F3(E,SQRT(MF2),1,UV,DV,US,DS,SS,FC,FB,GL)
       ELSEIF (IPDF.EQ.2) THEN
          SS  = CTQ5PDF(-3, E, SQRT(MF2))*E
          US  = CTQ5PDF(-1, E, SQRT(MF2))*E
          DS  = CTQ5PDF(-2, E, SQRT(MF2))*E
          GL  = CTQ5PDF( 0, E, SQRT(MF2))*E
          DV  = CTQ5PDF( 2, E, SQRT(MF2))*E
          UV  = CTQ5PDF( 1, E, SQRT(MF2))*E
c          write(6,*)'GL',GL
c          SV   = CTQ5PDF( 3, E, SQRT(MF2))*E
       ENDIF

* alpha_s using renormalisation scale
       AAS = ALPHA_S(sqrt(MR2)) / (4.D00*PI)
       AAS2 = AAS*AAS

C ... Note that UV = u - \bar{u} !
cpt think this is for GRV only
cpt       SF = UV + DV + 2.D00*US + 2.D00*DS + 2.D00*SS
cpt       SF2 = UV + 2.D00*US + 0.25D00 * (DV + 2.D00*DS + 2.D00*SS)
       SF = UV + DV + US + DS + 2.D00*SS
       SF2 = UV + US + 0.25D00 * (DV + DS + 2.D00*SS)
*
* ..The coefficient functions [expansion parameter alpha_s/(4 pi)]
* 
cpt check what is used?
       LOGFM = LOG (MF2/MH2)      
*
       Z  = X/E
       XI = Q2/MH2
       ETA = Q2 * (1.D0 - Z) / (4.D0 * Z * MH2) - 1.D0
*
C..    LO contribution:
       F2LO = AAS * XI / PI * GL * ( CTG0 (ETA,XI) + CLG0 (ETA,XI) )
       F2HI = F2LO

C..    The three NLO contributions:
       F2NLO1 = CTG1 (ETA,XI) + CLG1 (ETA,XI)
      1       + LOGFM * ( CTG1L (ETA,XI) + CLG1L (ETA,XI) )
       F2NLO1 = GL * F2NLO1  
       F2NLO2 = CTQ1 (ETA,XI) + CLQ1 (ETA,XI)
      1       + LOGFM * ( CTQ1L (ETA,XI) + CLQ1L (ETA,XI) )
       F2NLO2 = SF * F2NLO2  
       F2NLO3 = SF2 * ( DTQ1 (ETA,XI) + DLQ1 (ETA,XI) )

       F2NLO = F2NLO1 + F2NLO2 + F2NLO3
       F2NLO = 16.D00 * PI * XI * F2NLO
          
C..    Sum up LO and NLO
       F2HI = F2LO + AAS2 * F2NLO

       F2HI = F2HI / E
*
       RETURN
       END

* =====================================================================
*           F  2  L
* =====================================================================
* ..The convolution integrand
*
       FUNCTION F2LI (E)
*
       IMPLICIT DOUBLE PRECISION (A - H, L - Z)
*
* ..Input common-block

       PARAMETER ( PI = 3.1415926536D+00 )
*      
       COMMON / PARF2C / X, Q2, MF2, MR2, MH2, EH2

       COMMON / USEPDF/ IPDF
       INTEGER IPDF
*
* ..The parton densities at x = E (GL = x*g, as usual)
*
c       ISET = 2
c       CALL GRV98PA (ISET, E, MF2, UV, DV, US, DS, SS, GL)
       IF (IPDF.EQ.1) THEN
          FB  = 0D0
          FC  = 0D0
          CALL MRST2004F3(E,SQRT(MF2),1,UV,DV,US,DS,SS,FC,FB,GL)
       ELSEIF (IPDF.EQ.2) THEN
          SS  = CTQ5PDF(-3, E, SQRT(MF2))*E
          US  = CTQ5PDF(-1, E, SQRT(MF2))*E
          DS  = CTQ5PDF(-2, E, SQRT(MF2))*E
          GL  = CTQ5PDF( 0, E, SQRT(MF2))*E
          DV  = CTQ5PDF( 2, E, SQRT(MF2))*E
          UV  = CTQ5PDF( 1, E, SQRT(MF2))*E
c          SV   = CTQ5PDF( 3, E, SQRT(MF2))*E
       ENDIF

c use renormalisation scale for alpha_s
       AAS = ALPHA_S(sqrt(MR2)) / (4.D00*PI)
       AAS2 = AAS*AAS

C ... Note that UV = u - \bar{u} !
cpt we think this is just for GRV
cpt       SF = UV + DV + 2.D00*US + 2.D00*DS + 2.D00*SS
cpt       SF2 = UV + 2.D00*US + 0.25D00 * (DV + 2.D00*DS + 2.D00*SS)
       SF = UV + DV + US + DS + 2.D00*SS
       SF2 = UV + US + 0.25D00 * (DV + DS + 2.D00*SS)

*
* ..The coefficient functions [expansion parameter alpha_s/(4 pi)]
*
       LOGFM = LOG (MF2/MH2)      
*
       Z  = X/E
       XI = Q2/MH2
       ETA = Q2 * (1.D0 - Z) / (4.D0 * Z * MH2) - 1.D0

C..    LO contribution:
       F2LO = AAS * XI / PI * GL * CLG0 (ETA,XI)
 
C..    The three NLO contributions:
       F2NLO1 = CLG1 (ETA,XI) + LOGFM * CLG1L (ETA,XI) 
       F2NLO1 = GL * F2NLO1  
       F2NLO2 = CLQ1 (ETA,XI) + LOGFM * CLQ1L (ETA,XI)
       F2NLO2 = SF * F2NLO2  
       F2NLO3 = SF2 * DLQ1 (ETA,XI)

       F2NLO = F2NLO1 + F2NLO2 + F2NLO3
       F2NLO = 16.D00 * PI * XI * F2NLO
       
C..    Sum up all   
       F2LI = F2LO + AAS2 * F2NLO
       F2LI = F2LI / E
*
       RETURN
       END
*
* =====================================================================
*
*
* ..A simple gluon density (Les Houches 2001 benchmark input)
*
cpt       SUBROUTINE GLUON (X, Q2, GL)
cpt       IMPLICIT DOUBLE PRECISION (A - Z)
cpt       INTEGER ISET
*
C       GL = 1.7D0 * X**(-0.1D0) * (1.D0-X)**5

cpt       ISET = 2
cpt       CALL GRV98PA (ISET, X, Q2, UV, DV, US, DS, SS, GL)

*
cpt       RETURN
cpt       END

* =====================================================================

      double precision function alpha_s(q)
*     
*     
      IMPLICIT DOUBLE PRECISION (A-H, L-Z)

      DATA PI/3.145926535D00/

      COMMON / USEPDF/ IPDF
      INTEGER IPDF
cpt must be 3-flavours and have correct lambda QCD
* GRV
cpt      DATA XLA/0.25D00/
*     

      IF (IPDF.EQ.1) THEN
         XLA=0.407D00
c         XLA=0.393D00
      ELSEIF (IPDF.EQ.2) THEN
         XLA=0.395D00
      ENDIF


cpt      NF=4.D00 
      NF=3.D00 
      B0=11.D00-2.D00*NF/3.D00
      B1=102.D00-38.D00*NF/3.D00
      SIL=LOG(Q**2/XLA**2)
      DOL=LOG(SIL)
*     
* nlo 
      alfsn=4.D00*PI*(1.D00-B1*DOL/(SIL*B0*B0))/(B0*SIL)
      alpha_s=alfsn
*     
      return
      end


C =====================================================================
*     Simple parton densities without evolution 
*     for testing purposes

      SUBROUTINE TOYPART(X, Q2, UV, DV, US, DS, SS, GL)
      IMPLICIT DOUBLE PRECISION (A - Z)

      UV = 5.1072D00 * x**(0.8D00) * (1.D00 - x)**(3.D00)
      DV = 3.06432D00 * x**(0.8D00) * (1.D00 - x)**(4.D00)
      GL = 1.7D00 * x**(-0.1D00) * (1.D00 - x)**(5.D00)
      DS = 0.1939875D00 * x**(-0.1D00) * (1.D00 - x)**(6.D00)
      US = (1.D00 - x) * DS
      SS = 0.2D00 * (US + DS)

      return
      end

*
* =====================================================================
C       include 'asymp.f' 
C       include 'cg0f2lh.f'
C       include 'cg1f2lh.f'
C       include 'cq1f2lh.f'
C       include 'daind.f'
C       include 'dilog.f'
C       include 'hblg_ca.f'
C       include 'hblq_h.f'
C       include 'hbtg_ca.f'
C       include 'hbtq_h.f'
C       include 'hlg_ca.f'
C       include 'hlg_cf.f'
C       include 'hlq_h.f'
C       include 'htg_ca.f'
CC       include 'htg_cf.f'
C       include 'htq_h.f'
C       include 'locate.f'
C       include 'thresh.f'
C       include 'hlq_l.f'
C       include 'htq_l.f'
