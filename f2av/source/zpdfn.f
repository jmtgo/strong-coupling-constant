C     -------------------------------------------------
      REAL FUNCTION ZPDFN(X,Q2,PDFG)
C     --------------------------------------------------

 
      REAL X,Q2,PDFG(161,161)
      REAL Q2G1,Q2G2,XG1,XG2

      COMMON/ZPDFNCMN/Q2G1(0:159),Q2G2(0:159),XG1(0:159),XG2(0:159)

      LOGICAL LFIRST
      DATA    LFIRST /.TRUE./
      SAVE    LFIRST
     
      IF (Q2.LT.0.39.OR.X.LT.1E-6.OR.X.GT.0.977) THEN
        ZPDFN=0.0
        RETURN
      ENDIF

C
      IQ2L = 0
      IQ2H = 0
      IXL  = 0
      IXH  = 0

      TQ   = 0.
      TX   = 0.

      PDFH  = 0.
      PDFL  = 0.
      ZPDFN = 0.

C     Calculate grids points
      IF (LFIRST) THEN

       DO I=0,159
        IF(I.LE.71)THEN
          Q2G1(I)=10**(5.*I/120.)
          Q2G2(I)=10**(5.*(I+1)/120.)
        ELSEIF(I.EQ.72)THEN
          Q2G1(I)=10**(5.*I/120.)
          Q2G2(I)=10**(2.*(I-71)/88. +3.)
        ELSE
          Q2G1(I)=10**(2.*(I-72)/88. +3.)
          Q2G2(I)=10**(2.*(I-71)/88. +3.)
        ENDIF
       ENDDO

       DO I=0,159
        IF(I.LE.79)THEN
          XG1(I)=10**(5.*I/133.-5.)
          XG2(I)=10**(5.*(I+1)/133.-5.)
        ELSEIF(I.EQ.80)THEN
          XG1(I)=10**(5.*I/133.-5.)
          XG2(I)=10**(2./80.*(I-79)-2.01)
        ELSE
          XG1(I)=10**(2./80.*(I-80)-2.01)
          XG2(I)=10**(2./80.*(I-79)-2.01)
        ENDIF
       ENDDO
       LFIRST=.FALSE.

      ENDIF

C     Bilinear interpolation 

      DO I=0,159
        IF (Q2.GE.Q2G1(I).AND.Q2.LT.Q2G2(I)) THEN
          IQ2L=I+1
          IQ2H=I+2
          TQ=(Q2-Q2G1(I))/(Q2G2(I)-Q2G1(I))
          GOTO 777
        ENDIF
      ENDDO

 777  CONTINUE

      DO I=0,159
        IF (X.GE.XG1(I).AND.X.LT.XG2(I)) THEN
          IXL=I+1
          IXH=I+2
          TX=(X-XG1(I))/(XG2(I)-XG1(I))
          GOTO 888
        ENDIF
      ENDDO
    
 888  CONTINUE 

      PDFH=(1-TX)*PDFG(IXL,IQ2H)+TX*PDFG(IXH,IQ2H)
      PDFL=(1-TX)*PDFG(IXL,IQ2L)+TX*PDFG(IXH,IQ2L)
      ZPDFN=(1-TQ)*PDFL+TQ*PDFH

      RETURN 
      END
C     --------------------------------------------------
