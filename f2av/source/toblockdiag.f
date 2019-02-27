      Subroutine ToBlockDiag(diag,last,corr,box,dfl,dxf3,diagst,dflst)
C---------------------------------------------------------------
C
C Created on 5 July 2004 by SG.
C
C Modified Nov 10 by SG to add DFL, DXF3 matrices
C
C Transform problem to block-diagonal matrix. Matrices "box" and "last"
C are modified
C
C---------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 dfl(NF2MAX),dxf3(NF2MAX),dflst(NF2MAX)
      real*8 Box(NSYSTMAX,NSYSTMAX) ! Original syst --> covariance matrix
      real*8 DiagST(NF2MAX)
C Local:
      real*8 CorrT(NSYSTMAX,NF2MAX)
      integer if2,isys,j,iexp,isys1,isys2,irea
      real*8    Coef,SUM,chi2,chi2loc,sum1,sum2,pull
      real*8 Box2(NSYSTMAX,NSYSTMAX) ! Covariance -> Variance -> eigenvectors
     $     ,WWW(NSYSTMAX)

C Check:
      real*8 Box3(NSYSTMAX,NSYSTMAX) ! Variance
     $     ,Box4(NSYSTMAX,NSYSTMAX)  ! eigenvectors^-1 (= trans)
     $     ,Box5(NSYSTMAX,NSYSTMAX)  ! Diagonal with eigenvalues of variance mat
     $     ,Box6(NSYSTMAX,NSYSTMAX)  ! Check 1
     $     ,Box7(NSYSTMAX,NSYSTMAX)  ! Check 2 - should be equal to Variance (box3)


C Rotated systematics:
      real*8 Corr2(NSYSTMAX,NF2MAX)

      real y
      integer ndf
      real*8 work(100000)
      integer iwork(100000)
      integer ifail
      character*132 ctmp
      
      real sterr,uncerr
      integer imeas
      real uncor_syst
C Check:
      real*8 chiave1,chiave2,chistd
      integer ick
      real*8 F2Ran(NF2MAX),SYSRAN(NSYSTMAX),SysRan2(NSYSTMAX)

      character*8 ctmpname

C---------------------------------------------------------------------
      

C
C === Find new Central Values and Uncertainties for the Systematics
C


C Copy corr to corrT:
      do if2=1,NSFSEC
         do isys=1,NSYSTOT
            CorrT(isys,if2) = Corr(isys,if2)
         enddo
      enddo

C      IDEBUG =10

C Pre-process for FL, xF3:

      if (IAVMETH.ne.0) then
         Call FLXF3PROC(diag,last,corr,corrT,dfl,dxf3,diagst,dflst)
      endif

C Zero entries under diagonal F2 entries:

      if (IDEBUG.gt.4) then
         print *,'in toblock'
         do j=1,nsystot
            print *,'systematic=',j,sysnam(j)
            print '(10F8.4)',(box(isys,j),isys=1,nsystot)
         enddo
         print *,(last(j),j=1,nsfsec+nsystot)
         print *,(diag(j),j=1,nsfsec)
         print *,'corrT:'
         do if2=1,nsfsec
            print *,(corrT(j,if2),j=1,nsfsec)
         enddo
         print *,'corr:'
         do if2=1,nsfsec
            print *,(corr(j,if2),j=1,nsfsec)
         enddo
      endif


      do if2=1,NSFSEC
         do isys=1,NSYSTOT
            Coef = - Corr(isys,if2)/diag(if2)
            if (IDEBUG.gt.4) then
               print *,if2,isys,coef
            endif
            if (Coef.ne.0) then
               do j=1,NSYSTOT
                  box(j,isys) = box(j,isys) + corrT(j,if2)*Coef
               enddo
               last(NSFSEC+isys) = last(NSFSEC+isys) + last(if2)*Coef

            endif
         enddo
      enddo

      if (IDEBUG.gt.4) then
         do j=1,nsystot
            print '(10F8.4)',(box(isys,j),isys=1,nsystot)
         enddo
         print *,(last(nsfsec+j),j=1,nsystot)
      endif
C
C Invert box matrix corresponding to systematic uncertainties:
C Here we obtain new systematic central values/uncertainties
C
      Call DEQINV(nsystot,box,nsystmax,work,ifail,1,last(nsfsec+1))      

      if (IDEBUG.gt.4) then
         print *,ifail
         do j=1,nsystot
            print '(10F8.4)',(box(isys,j),isys=1,nsystot)
         enddo
         print *,(last(nsfsec+j),j=1,nsystot)
      endif

      if (NSTATADJ.eq.0 .or. IDEBUG.gt.2) then
         print *,' '
         print *,'Fitted systematics:'
         print *,' '
         do j=1,nsystot
            errsyst(j) = sqrt(box(j,j))
            ctmp = sysnam(j)
            ctmp = ctmp(1:index(ctmp,' '))
            print '(i3,''  '',A32,2F10.4)',j,ctmp,last(nsfsec+j)
     $           ,errsyst(j)
         enddo
      endif

C 15 July SG: Store the systematics:

      do j=1,NSYSTOT
         SYSSH(j) = last(nsfsec+j)
      enddo

C
C ==== CALCULATE AVERAGE X-SECTIONS
C

C
C Now we have systematics (with errors !), get F2s:
C
      do if2=1,NSFSEC
         SUM = Last(if2)
         do isys=1,NSYSTOT
            Sum = Sum - SYSSH(isys)*CorrT(isys,if2)
         enddo
         F2VAVE(if2) = SUM/DIAG(if2)

      enddo



C Check if we want to do the following or just bail out:
      if (NSTATADJ.eq.0 .or. IDEBUG.gt.2) then
      else
C Do nothing otherwise 
         Return
      endif


C 11 Jan 06: Zero close to zero correlations to avoid strange results
C            for a single dataset
      do isys1=1,nsystot
         do isys2=1,nsystot
C            print *,box(isys1,isys2)
            if (abs(box(isys1,isys2)).lt.1.E-10) then
               box(isys1,isys2) = 0
            endif
         enddo
         if (abs(box(isys1,isys1)-1.0).lt.1.E-10) then
            box(isys1,isys1) = 1.0D0
         endif
      enddo

C 8 july 04: Save  "box" in box2:

      do isys1=1,nsystot
         do isys2=1,nsystot
            box2(isys1,isys2) = box(isys1,isys2)
C     $           /sqrt(box(isys1,isys1)*box(isys2,isys2))
         enddo
      enddo

C Invert box2 to get variance matrix:      
      Call DINV(nsystot,box2,nsystmax,work,ifail)

C Save:
      do isys1=1,nsystot
         do isys2=1,nsystot
           box3(isys1,isys2) = box2(isys1,isys2) 
         enddo
      enddo      


C
C ====      DIAGONALIZE SYSTEMATICS ===  
C

C Get eigenvectors and eigenvalues:

      Call DSYEVD('V','U',nsystot,Box2,nsystmax,WWW,
     $     work,100000,iwork,100000,ifail)


      if (mod(iout/1000,10).eq.1) then
C Print the eigenvalues and eigenvectors
         open (51,file='output/eigvalues.dat',status='unknown')
         do isys1=1,nsystot
            write (51,'(E16.8)') WWW(isys1)
         enddo         
         close (51)
         open (51,file='output/eigvectors.dat',status='unknown')
         do isys1=1,nsystot
            write (51,'(100E16.8)') (box2(isys2,isys1),isys2=1,nsystot)
         enddo         
         close (51)
      endif

      if (ifail.ne.0) then
         print *,
     $   'Failed to find eigenvectors for systematic sourses ifail='
     $        , ifail
         STOP
      endif


C Store in Box4:
      do isys1=1,nsystot
         do isys2=1,nsystot
           box4(isys1,isys2) = box2(isys1,isys2) 
         enddo
      enddo      

C Find inverse:
      Call DINV(nsystot,box4,nsystmax,work,ifail)



      if (IDEBUG.ge.3) then
         print *,'P:'
         do isys1=1,nsystot
            print '(100F6.2)',(box2(isys1,isys2),isys2=1,nsystot)
         enddo

         print *,'P-1:'
         do isys1=1,nsystot
            print '(100F6.2)',(box4(isys1,isys2),isys2=1,nsystot)
         enddo
      endif

C Diagonal matrix:
      do isys=1,nsystot
         box5(isys,isys) = WWW(isys)
      enddo


C Build P*D*P^-1 product
      do isys1=1,nsystmax 
         do isys2=1,nsystmax
            do isys=1,nsystmax
               box6(isys1,isys2) = box6(isys1,isys2) +
     $              box2(isys1,isys)*box5(isys,isys2)
            enddo
         enddo
      enddo
      
      do isys1=1,nsystmax 
         do isys2=1,nsystmax
            do isys=1,nsystmax
               box7(isys1,isys2) = box7(isys1,isys2) +
     $              box6(isys1,isys)*box4(isys,isys2)
            enddo
         enddo
      enddo



      if (IDEBUG.ge.1) then

         print *,'Systematics correlation matrix:'
         do isys1=1,nsystot
            print '(100F6.2)',(box(isys1,isys2)
     $  /sqrt(box(isys1,isys1)*box(isys2,isys2)),isys2=1,nsystot)
         enddo

         print *,'(Reduced) systematics variance matrix:'
         do isys1=1,nsystot
            print '(100F6.2)',(box3(isys1,isys2),isys2=1,nsystot)
         enddo
      endif

      if (IDEBUG.ge.3) then
         print *,'Check:'
         do isys1=1,nsystot
            print '(100F6.2)',(box7(isys1,isys2),isys2=1,nsystot)
         enddo
      endif

C
C Get rotated systematic matrix:
C
      do if2=1,NSFSEC
         do isys1=1,NSYSTOT
            do isys=1,NSYSTOT
               Corr2(isys1,if2) = Corr2(isys1,if2) +
     $              CorrT(isys,if2)*Box2(isys,isys1)
            enddo
C            print *,if2,isys1,corr2(isys1,if2),corr(isys1,if2)
         enddo
      enddo



C
C Get NEW syst. errors for F2:
C
      do if2=1,NSFSEC
         SUM = 0
         do isys1=1,NSYSTOT
            do isys2=1,NSYSTOT
               SUM = SUM + CorrT(isys1,if2)*
     $              CorrT(isys2,if2)/diag(if2)**2
     $              * box(isys1,isys2)
            enddo
         enddo
         F2EsyAve(if2) = sqrt(Sum)
         if (IDEBUG.gt.3) then
            print *,if2,sqrt(Sum)
         endif
      enddo

C CHECK:

      if (IDEBUG.gt.3) then
         do if2=1,NSFSEC
            SUM = 0.
            do isys=1,NSYSTOT
               Sum = Sum + 
     $              Corr2(isys,if2)**2/diag(if2)**2/box5(isys,isys)
            enddo
            print *,if2,sqrt(Sum)
         enddo
      endif

C
C Get stat. and total errors for F2:
C

      do if2=1,NSFSEC
      enddo
      

      do if2=1,NSFSEC
C proper way:
         STERR = 0.
         UNCERR = 0.
C Stat. and uncorr errors:
         F2EstAve(if2)     = sqrt(1/diag(if2))
C True stat. errors:
         do imeas =1,NMEAS
            if (F2TAB(if2,imeas).ne.0) then
               STERR = STERR  +
     $              F2ETAB_STA(if2,imeas)**2/F2ETAB(if2,imeas)**4
               UNCERR = UNCERR +
     $              (F2ETAB(if2,imeas)**2-F2ETAB_STA(if2,imeas)**2)
     $              /F2ETAB(if2,imeas)**4
            endif
         enddo
         STERR = sqrt(STERR/(DIAG(if2)**2))
         UNCERR = sqrt(UNCERR/(DIAG(if2)**2))
         if (IAVMETH.ne.0) then
            F2EstAveTrue(if2) = sqrt(1/diagst(if2)) 
         else
            F2EstAveTrue(if2) = STERR ! sqrt(1/diagst(if2)) 
         endif
         F2EAVE(if2)       = sqrt(F2EsyAve(if2)**2+F2EstAve(if2)**2)         
      enddo


C
C == GET CHI2/DOF, prepare for output
C


      if (IAVMETH.eq.0) then
C Time for chi2:

         open(51,file='output/chi2map.dat',status='unknown')
         chi2 = 0.0
         ndf  =  0
         do iexp=1,NMEAS
            do if2=1,NF2TOT
               if (F2TAB(if2,iexp).ne.0) then
                  ndf = ndf + 1
                  sum = F2TAB(if2,iexp)
                  do isys=1,NSYSTOT
                     sum = sum + SYSTAB(isys,if2,iexp)*SYSSH(isys)
                  enddo
                  chi2loc = (F2VAVE(if2)-sum)
     $                 /F2ETAB(if2,iexp)
                  chi2 = chi2 + chi2loc**2
C For pulls, divide by err^2 - err_ave^2 instead of err^2:

                  if (NMEASF2(if2).gt.1) then
                     pull = (F2VAVE(if2)-sum)/
     $                  sqrt(abs(F2ETAB(if2,iexp)**2-F2EstAve(if2)**2))
                  else
                     pull = 0
                  endif
                  
                  write (51,'(i4,i4,5F14.6)') 
     $                 iexp,nmeasf2(if2),Q2TAB(if2),XTAB(if2),
     $                 pull,f2tab(if2,iexp),Sum
               endif
            enddo
         enddo

         
         close (51)




C
C INPUT DATA <==> AVERAGE CONSISTENCY CHECK
C

C 9 Aug 04: Check that the average is equivalent to the input data.
C For that generate ICHECK random F2 and Systematic variations 
C
         if (ICHECK.ne.0) then
            
            print *,' '
            print *,'Check Chi2 for several points'
            print *,' '
            print '(A4,3A16,2A16)',' Check ',' Std Chi2 ',
     $        'Ave1 Chi2 ','Ave2 Chi2 ',
     $           'Ave1/Std-1 ','Ave2/Std-1 '
         

            do ick = 0,ICHECK
C
C Generate random F2 and Systematic:
C     
               Call F2RND(F2ran,Sysran,ick)

C Rotate the systematic matrix:
               do isys1=1,NSYSTOT
                  SysRan2(isys1) = 0
                  do isys2=1,NSYSTOT
                     SysRan2(isys1) = SysRan2(isys1) +
     $                    Box2(isys2,isys1)*
     $                    (SysRan(isys2)-SysSh(isys2))
                  enddo
                  SysRan2(isys1) = SysRan2(isys1)
               enddo

C
C Zero Chi2:
C     
               chistd = 0.
               chiave1 = chi2
               chiave2 = chi2

C Sum F2 contribution:
               do if2=1,NF2TOT
C     Get chi2 for this point using standard and average tables:
                  do iexp=1,NMEAS
                     if (F2TAB(if2,iexp).ne.0) then
                        sum = F2TAB(if2,iexp)
                        do isys=1,NSYSTOT
                           sum = sum 
     $                          + SYSTAB(isys,if2,iexp)*SYSRan(isys)
                        enddo
                        chistd = chistd + 
     $                       (Sum-F2Ran(if2))**2
     $                       /F2ETAB(if2,iexp)**2
                     endif
                  enddo

C Get for average (2 different syst. representations):
                  sum = F2VAVE(if2)              
                  sum1 = 0
                  sum2 = 0

                  do isys=1,NSYSTOT
                     sum1 = sum1 
     $                    - CORR(isys,if2)/DIAG(if2)
     $                    *(SYSRan(isys)-SYSSH(isys))
                     sum2 = sum2 
     $                    - CORR2(isys,if2)/DIAG(if2)
     $                    *SYSRan2(isys)
                  enddo

                  
                  chiave1 = chiave1 + (sum + sum1-F2Ran(if2))**2
     $                 /F2EstAve(if2)**2
                  
                  chiave2 = chiave2 + (sum + sum2-F2Ran(if2))**2
     $                 /F2EstAve(if2)**2

               enddo

C Sum systematics contribution:
               do isys=1,NSYSTOT
                  chistd = chistd + SYSRan(isys)**2
               enddo


               do isys1=1,NSYSTOT
                  do isys2=1,NSYSTOT
                     chiave1 = chiave1 + (SYSRan(isys1)-SYSSH(isys1))*
     $                    (SYSRan(isys2)-SYSSH(isys2))
     $                    *Box3(isys1,isys2)
                  enddo
               enddo


               do isys=1,NSYSTOT
                  chiave2 = chiave2 + SysRan2(isys)**2*Box5(isys,isys)
                  
               enddo
            
               print '(I4,3E16.6,2E16.6)'
     $              ,ick,chistd,chiave1,chiave2,chiave1/chistd-1.
     $              ,chiave2/chistd-1.

            enddo
         

         endif



         if (LOGAVE) then
C
C Back from LOG(S_red) to S_red:
C
            do if2=1,NF2TOT
               f2vave(if2) = exp(f2vave(if2))
               f2estave(if2) = f2estave(if2)*f2vave(if2)
               F2EstAveTrue(if2) =   F2EstAveTrue(if2) *f2vave(if2) 
               f2esyave(if2) = f2esyave(if2)*f2vave(if2)
               f2eave(if2)   = f2eave(if2)  *f2vave(if2)
               do isys=1,NSYSTOT
                  CORR2(isys,if2) = CORR2(isys,if2)*f2vave(if2)
               enddo
            enddo
         endif

         print *,' '
         print *,'Output F2:'
         print *,' '
         
         do irea=1,NREACT
            print *,' '
            print '('' Reaction = '',i5)' ,RLIST(irea)
            print *,' '
            print '(7A10)',' Q2  ',' X  ','  Y  ',' X-sect',' E-Uncor',
     $           ' E-Corr',' E-Total'
            if (rlist(irea).lt.10) then
               write (ctmpname,'(i1)') rlist(irea)
            elseif (rlist(irea).lt.100) then
               write (ctmpname,'(i2)') rlist(irea)
            elseif (rlist(irea).lt.1000) then
               write (ctmpname,'(i3)') rlist(irea)
            elseif (rlist(irea).lt.10000) then
               write (ctmpname,'(i4)') rlist(irea)
            elseif (rlist(irea).lt.100000) then
               write (ctmpname,'(i5)') rlist(irea)
            endif

            open (51,file='output/xsec'//ctmpname(1:len(ctmpname))
     $           ,status='unknown')
            do if2=1,NF2TOT
               if (REACTAB(if2).eq.RLIST(irea)) then
                  print '(F10.2,E14.4,5F10.4)',
     $                 q2tab(if2),xtab(if2),ytab(if2),f2vave(if2)
     $                 ,f2estave(if2),f2esyave(if2),f2eave(if2)
                  write (51, '(F10.2,E14.4,5F14.8)')
     $                 q2tab(if2),xtab(if2),ytab(if2),f2vave(if2)
     $                 ,f2estave(if2),f2esyave(if2),f2eave(if2)
               endif
            enddo
            close (51)
         enddo

         print *,' '
         print '(''Chi2 EXCLUDING systematics='',F10.4)'
     $        ,chi2

         do isys=1,NSYSTOT
            chi2 = chi2 + SYSSH(isys)**2
         enddo

         print *,' '
         print '(''TOTAL Chi2/ndf='',F10.4,''/'',i4)',chi2,ndf-nf2tot
         print *,' '

      else if (IAVMETH.eq.1 .or. IAVMETH.eq.10) then
C F2/FL output
         Call WriteF2FL(Corr2,Box,Box5,diag,diagst)
      endif

C
C ============        DATA OUTPUT     ==================
C
      if (IAVMETH.eq.0) then

C
C Also write out data in a "standard format" if requested:
C
         if (mod(IOUT/100,10).ne.0) then
            open (51,file='output/tab.dat',status='unknown')
            
            write 
     $      (51,'(''*  Total number of systematic sources *'')')
            write (51,'(I5)') NSYSTOT
            
C For options 1 and 2 write out the systematics correlation matrix:
            if (mod(IOUT/100,10).ne.3) then
               write (51,'(''* '')')
               write (51,
     $         '(''* Correlation matrix for systematic sources'')')
               write (51,'(''* '')')
               do isys1=1,nsystot
                  write (51,'(i4,100F10.6)') isys1,
     $                 (Box(isys1,isys2)
     $                 /sqrt(Box(isys1,isys1)*Box(isys2,isys2))
     $                 ,isys2=1,NSYSTOT)
               enddo
            endif

            if (mod(IOUT/100,10).eq.1) then
C
C Write in human readable output:
C
               write (51,'(''* '')')
               write (51,
     $ '(''*  All F2 and correlations (all errors are in percent) *'')')
               write (51,'(''* '')')
               write (51,'(''* Total number of data values: '')')
               write (51,'(''* '')')
               write (51,'(I5)') NF2TOT
               write (51,'(''* '')')
               write (51,
     $              '(''* Number of data points by reaction type '')')
               write (51,'(''* '')')
               write (51,'(100(I5,'':'',I5))') 
     $              (RLIST(irea),RCOUNT(irea),irea=1,NREACT)
               do irea=1,NREACT
                  write (51,'(''* Reaction type'')')
                  write (51,'(I5)') RLIST(irea)
                  write (51,'(''*'',A9,3A10,4X,100A10)')
     $                 ' Q2 ', ' X ',' Val ',' Uncor ',
     $                 (SYSNAM(isys),isys=1,NSYSTOT)
                  do if2=1,NF2TOT
                     if (REACTAB(if2).eq.RLIST(irea)) then
                        write (51,'(4F12.6,100F6.2)')
     $                       q2tab(if2),xtab(if2),f2vave(if2)
     $                       ,100.*f2estave(if2)/F2vave(if2),(
     $                       -100.*CORR(ISYSR(isys,irea),if2)/DIAG(if2)
     $                       *ErrSyst(ISYSR(isys,irea))
     $                       /f2vave(if2)
     $                       ,isys=1,NSYSR(irea))
                     endif
                  enddo
               enddo
            else if (mod(IOUT/100,10).eq.2) then
C
C Write in output for QCDFIT program, using original syst. sources:
C
               do irea=1,NREACT
 1717             format (i5,a8,2i3,i5,F10.3,1X,a1,100i2)
                  write (51,1717) RCOUNT(irea)
     $                 ,'AVE',
C     July 22 04: Write only relevant systematics:
     $                 NSYSR(irea),NSYSR(irea),RLIST(irea),
     $                 COMEXP,'O',(ISYSR(isys,irea),isys=1,NSYSR(irea))
                  do if2=1,NF2TOT
                     if (REACTAB(if2).eq.RLIST(irea)) then
 1718                   format (8E14.7)
                        
C For now Y uses CM for the first experiment:
                        y = ytab(if2) ! q2tab(if2)/(COMEXP*XTAB(if2))
                        uncor_syst = sqrt(f2estave(if2)**2
     $                       -f2estaveTrue(if2)**2)

                        write (51,1718) 
     $                       xtab(if2)
     $                       ,y
     $                       ,q2tab(if2)
     $                       ,f2vave(if2)
     $                       ,f2estaveTrue(if2)
     $                       ,uncor_syst
     $                       ,F2EAVE(if2)
 1719                   format (200F8.4)
c                        write (51,1719) (
c     $                       -100.*CORR(ISYSR(isys,irea),if2)/DIAG(if2)
c     $                       *ErrSyst(ISYSR(isys,irea))
c     $                       /f2vave(if2)
c     $                       ,isys=1,NSYSR(irea))
                        write (51,1719) (
     $                       -100.*CORR(isys,if2)/DIAG(if2)
     $                       *ErrSyst(isys)
     $                       /f2vave(if2)
     $                       ,isys=1,NSYSTOT)
                     endif
                  enddo
               enddo
            else if (mod(IOUT/100,10).eq.3) then
C     
C  Write in output for QCDFIT program, using orthogonal syst. sources:
C
               do irea=1,NREACT
                  write (51,1717) RCOUNT(irea)
     $                 ,'AVE',
     $                 NSYSTOT,NSYSTOT,RLIST(irea),
     $                 COMEXP,'O',(isys,isys=1,min(99,NSYSTOT))
                  
                  do if2=1,NF2TOT
                     if (REACTAB(if2).eq.RLIST(irea)) then
C For now Y uses CM for the first experiment: 
                        y = ytab(if2) ! q2tab(if2)/(COMEXP*XTAB(if2))
                        uncor_syst = sqrt(f2estave(if2)**2
     $                       -f2estaveTrue(if2)**2)

                        write (51,1718) xtab(if2),y,q2tab(if2),
     $                       f2vave(if2), 
     $                       F2estaveTrue(if2),
     $                       uncor_syst,
     $                       F2EAVE(if2)
                        write (51,1719) (
     $                       -100.*CORR2(isys,if2)/DIAG(if2)
     $                       /sqrt(Box5(isys,isys))
     $                       /f2vave(if2)
     $                       ,isys=1,NSYSTOT)
                        
                     endif
                  enddo               
               enddo            
            endif
         endif
      else
C
C F2 out
C
         open(51,file='output/f2.dat',status='unknown')
         do irea=1,NREACT
            write (51,1717) RCOUNT(irea)
     $           ,'F2',
     $           NSYSTOT,NSYSTOT,RLIST(irea),
     $           COMEXP,'O',(isys,isys=1,min(99,NSYSTOT))
            
            do if2=1,NF2TOT
               if (REACTAB(if2).eq.RLIST(irea)) then
C     For now Y uses CM for the first experiment: 
                  y = ytabmax(if2) ! q2tab(if2)/(COMEXP*XTAB(if2))
                  uncor_syst = sqrt(f2estave(if2)**2
     $                 -f2estaveTrue(if2)**2)
                  
                  write (51,1718) xtab(if2),y,q2tab(if2),
     $                 f2vave(if2), 
     $                 F2estaveTrue(if2),
     $                 uncor_syst,
     $                 F2EAVE(if2)
                  write (51,1719) (
     $                 -CORR2(isys,if2)/DIAG(if2)
     $                 /sqrt(Box5(isys,isys))
     $                 ,isys=1,NSYSTOT)
                        
               endif
            enddo               
         enddo            

         close (51)
C
C FL out
C
         open(51,file='output/fl.dat',status='unknown')

         do irea=1,NREACT
            write (51,1717) RCOUNT(irea)
     $           ,'FL',
     $           NSYSTOT,NSYSTOT,RLIST(irea),
     $           COMEXP,'O',(isys,isys=1,min(99,NSYSTOT))
            
            do if2=NF2TOT+1,2*NF2TOT
               if (REACTAB(if2-NF2TOT).eq.RLIST(irea)) then
C     For now Y uses CM for the first experiment: 
                  y = ytabmax(if2-NF2TOT) ! q2tab(if2)/(COMEXP*XTAB(if2))
                  uncor_syst = sqrt(f2estave(if2)**2
     $                 -f2estaveTrue(if2)**2)
                  
                  write (51,1718) xtab(if2-NF2TOT),y,q2tab(if2-NF2TOT),
     $                 f2vave(if2), 
     $                 F2estaveTrue(if2),
     $                 uncor_syst,
     $                 F2EAVE(if2)
                  write (51,1719) (
     $                 -CORR2(isys,if2)/DIAG(if2)
     $                 /sqrt(Box5(isys,isys))
     $                 ,isys=1,NSYSTOT)
                        
               endif
            enddo               
         enddo            
         
         close (51)
         
      endif
C---------------------------------------------------------------------
      end
