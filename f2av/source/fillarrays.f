      Subroutine FillArrays(diag,last,corr,box,dfl,dxf3,diagst,dflst)
C---------------------------------------------------------------
C
C Created on 5 July 2004 by SG.
C
C Modified 10 Nov 2007 by SG: add extra dfl and dxf3 arrays
C Modified  5 Dec 2008 by SG: add diagst for separated pure stat. errors
C
C Fill axillary arrays for the matrix inversion
C
C---------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 dfl(NF2MAX),dxf3(NF2MAX),dflst(NF2MAX)
      real*8 Box(NSystmax,NSystmax)
      real*8 diagst(NF2MAX)
C Local:
      integer isys,if2,iexp,j,i,isys1,isys2
      integer noff
      real*8 coef, erro
C---------------------------------------------------------------------

C Count number of SF / X-section points to fit 

      NSFSEC = NF2TOT
      if (MOD(IAVMETH,10).eq.1) then
C FL extraction
         NSFSEC = NSFSEC + NF2TOT   ! add FL points
      endif
      if (MOD(IAVMETH,100)/10.eq.1) then
C xF3 extraction
         NSFSEC = NSFSEC + NF2TOT   ! add xF3 points
      endif

      print *,'Total number of SF/X-section points=',NSFSEC,NF2TOT


C Zero everything:
      do if2=1,NSFSEC
         diag(if2)   = 0.0         
         diagst(if2) = 0.0
         do isys=1,NSYSTOT
            corr(isys,if2) = 0.0
         enddo
      enddo

      do if2=1,NF2TOT
         dfl(if2)  = 0.
         dxf3(if2) = 0.
         dflst(if2) = 0.
      enddo

      do if2=1,NSFSEC+NSYSTOT
         last(if2) = 0.0
      enddo

      do isys1=1,NSYSTOT
         do isys2=1,NSYSTOT
            box(isys1,isys2) = 0.0
         enddo
      enddo


C
C Start XS/F2 part:
C

C Diagonal XS/F2:
      do if2=1,NF2TOT
         do iexp=1,NMEAS
            if (F2TAB(if2,iexp).ne.0) then
               diag  (if2) = diag  (if2) + 1/F2ETAB(if2,iexp)**2
               diagst(if2) = diagst(if2) + 1/F2ETAB_STA(if2,iexp)**2
            endif
         enddo
      enddo

C Last column:
      do if2=1,NF2TOT
         do iexp=1,NMEAS
            if (F2TAB(if2,iexp).ne.0) then
               last(if2) = last(if2) + 
     $              F2TAB(if2,iexp)/F2ETAB(if2,iexp)**2
            endif
         enddo
      enddo


C Correlation:
      do if2=1,NF2TOT
         do iexp=1,NMEAS
            do isys=1,NSYSTOT
               if (F2TAB(if2,iexp).ne.0) then
                  CORR(isys,if2) = CORR(isys,if2) 
     $                 -SYSTAB(isys,if2,iexp)
     $                 /F2ETAB(if2,iexp)**2
               endif
            enddo
         enddo
      enddo

      noff = NF2TOT
C
C FL part:
C      
      if (IAVMETH.eq.1) then
C Diagonal FL:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  diag(if2+noff) = diag(if2+noff)+
     $                 YFACTFL(if2,iexp)**2/F2ETAB(if2,iexp)**2
                  diagst(if2+noff) = diagst(if2+noff)+
     $                 YFACTFL(if2,iexp)**2/F2ETAB_STA(if2,iexp)**2
C                  print *,'bb',diag(if2+noff) ,diagst(if2+noff)

C                  print *,if2,iexp,YFACTFL(if2,iexp)
               endif
            enddo
C            print *,'diag fl',if2,diag(if2+noff)
         enddo

C Last column:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  last(if2+noff) = last(if2+noff) -
     $    F2TAB(if2,iexp)*YFACTFL(if2,iexp)/F2ETAB(if2,iexp)**2
               endif
            enddo
C            print *,'last fl',if2,last(if2+noff)
         enddo

C Correlation:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               do isys=1,NSYSTOT
                  if (F2TAB(if2,iexp).ne.0) then
                     CORR(isys,if2+noff) = CORR(isys,if2+noff) 
     $                    +SYSTAB(isys,if2,iexp)*YFACTFL(if2,iexp)
     $                    /F2ETAB(if2,iexp)**2
C                     print *,'h',F2ETAB(if2,iexp),SYSTAB(isys,if2,iexp)
                  endif
               enddo
            enddo
C            print *,'corr fl',(corr(isys,if2+noff),isys=1,nsystot)
         enddo

C
C F2-FL diagonal:
C
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  dfl(if2) = dfl(if2) -
     $                 YFACTFL(if2,iexp)/F2ETAB(if2,iexp)**2
                  dflst(if2) = dflst(if2) -
     $                 YFACTFL(if2,iexp)/F2ETAB(if2,iexp)**2
               endif
            enddo
C            print *,'dfl',if2,dfl(if2)
         enddo

         if (LFLLIMITS) then
C add extra FL constraint to mimic  0 < FL < F2 in a "soft" way
            coef = RAVERAGE/(RAVERAGE + 1)
            erro = RERROR/(RERROR + 1)
            do if2=1,nf2tot
C F2 diag:
               diag(if2) = diag(if2) + coef**2/erro**2
               diagst(if2) = diagst(if2) + coef**2/erro**2
C FL diag:
               diag(if2+noff) = diag(if2+noff) + 1/erro**2
               diagst(if2+noff) = diagst(if2+noff) + 1/erro**2
C F2-FL diag:
               dfl(if2) = dfl(if2) - coef/erro**2
               dflst(if2) = dflst(if2) - coef/erro**2
            enddo
         endif
         
         noff = noff + NF2TOT
      endif

C
C xF3 part:
C
      if (IAVMETH.eq.10) then


C Diagonal xF3:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  diag(if2+noff) = diag(if2+noff)+
     $                 YFACTF3(if2,iexp)**2/F2ETAB(if2,iexp)**2
                  diagst(if2+noff) = diagst(if2+noff)+
     $                 YFACTF3(if2,iexp)**2/F2ETAB_STA(if2,iexp)**2
C                  print *,'bb',diag(if2+noff) ,diagst(if2+noff)

C                  print *,if2,iexp,YFACTF3(if2,iexp)
               endif
            enddo
c            print *,'diag xf3',if2,diag(if2+noff)
         enddo

C Last column:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  last(if2+noff) = last(if2+noff) -
     $    F2TAB(if2,iexp)*YFACTF3(if2,iexp)/F2ETAB(if2,iexp)**2
               endif
            enddo
c            print *,'last xf3',if2,last(if2+noff)
         enddo

C Correlation:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               do isys=1,NSYSTOT
                  if (F2TAB(if2,iexp).ne.0) then
                     CORR(isys,if2+noff) = CORR(isys,if2+noff) 
     $                    +SYSTAB(isys,if2,iexp)*YFACTF3(if2,iexp)
     $                    /F2ETAB(if2,iexp)**2
C                     print *,'h',F2ETAB(if2,iexp),SYSTAB(isys,if2,iexp)
                  endif
               enddo
            enddo
C            print *,'corr xf3',(corr(isys,if2+noff),isys=1,nsystot)
         enddo

C
C F2-xF3 diagonal:
C
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               if (F2TAB(if2,iexp).ne.0) then
                  dfl(if2) = dfl(if2) -
     $                 YFACTF3(if2,iexp)/F2ETAB(if2,iexp)**2
                  dflst(if2) = dflst(if2) -
     $                 YFACTF3(if2,iexp)/F2ETAB(if2,iexp)**2
               endif
            enddo
C            print *,'dfl',if2,dfl(if2)
         enddo


         if (LFLLIMITS) then
C add extra xF3 constraint to mimic  0 < xF3 < F2 in a "soft" way
            coef = RAVERAGE/(RAVERAGE + 1)
            erro = RERROR
            do if2=1,nf2tot
C F2 diag:
               diag(if2) = diag(if2) + coef**2/erro**2
               diagst(if2) = diagst(if2) + coef**2/erro**2
C xF3 diag:
               diag(if2+noff) = diag(if2+noff) + 1/erro**2
               diagst(if2+noff) = diagst(if2+noff) + 1/erro**2
C F2-FL diag:
               dfl(if2) = dfl(if2) - coef/erro**2
               dflst(if2) = dflst(if2) - coef/erro**2
            enddo
         endif

         noff = noff + NF2TOT

      endif



C last column correlation:
      do iexp=1,NMEAS
         do isys=1,NSYSTOT
            do if2=1,NF2TOT 
               if (F2TAB(if2,iexp).ne.0) then
                  last(Noff+isys) = last(Noff+isys) 
     $                 - SYSTAB(isys,if2,iexp)
     $                 * F2TAB(if2,iexp)/F2ETAB(if2,iexp)**2
               endif
            enddo
         enddo
      enddo



C Box:

C 5 July 2004 SG: for now assume no correlation between experiments:      

      do iexp=1,NMEAS
         do isys1=1,NSYSTOT
            do isys2=1,NSYSTOT
               do if2=1,nf2tot
                  if (F2TAB(if2,iexp).ne.0) then
                     box(isys1,isys2) = box(isys1,isys2)
     $                    + SYSTAB(isys1,if2,iexp)
     $                    * SYSTAB(isys2,if2,iexp)
     $                    / F2ETAB(if2,iexp)**2
                  endif
               enddo
            enddo
         enddo
      enddo


C Add 1 to box diagonal:

      do j=1,NSYSTOT
         BOX(j,j) = BOX(j,j) + 1.0
      enddo


      if (IDEBUG.gt.4) then
         print *,'cor'
         do if2=1,NF2tot
            print '(20F10.4)',(corr(j,if2),j=1,NSYSTOT)
         enddo
         print *,'box'
      
         do j=1,NSYSTOT
            print '(20F10.4)',(box(i,j),i=1,NSYSTOT)
         enddo
         print *,'last'
         print '(6F12.4)',(last(j),j=1,NF2TOT+NSYSTOT)
      endif
C---------------------------------------------------------------------
      end
