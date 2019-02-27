      Subroutine WriteF2FL(corr2,box,box5,diag,diagst)
C------------------------------------------------------
C
C  Write output for F2, FL and/or xF3 fit
C
C
C------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      real*8 corr2(NSYSTMAX,NF2MAX),box(NSYSTMAX,NSYSTMAX)
     $     ,box5(NSYSTMAX,NSYSTMAX)
     $     ,diag(NF2MAX),diagst(NF2MAX),chi2,sred
      integer if2,nsf,ndf,iexp,isys,isys1,isys2
      real*8 sum,f2,fl,chi2loc
C---------------------------------------------------------
      nsf = 1
      print *,' '
      print *,'Correlation of systematics'
 1    format (A4,50(I4))
 2    format (I4,50(I4))
      print 1,' ',(isys,isys=1,nsystot)
      do isys1=1,nsystot
         print 2,isys1,(int(100.*box(isys1,isys2)
     $        /sqrt(box(isys1,isys1)*box(isys2,isys2))),isys2=1,nsystot)
      enddo
      print *,' '
      print *,' '
      if (mod(IAVMETH,10).eq.1) nsf = nsf+1
      if (mod(IAVMETH,100)/10.eq.1) nsf = nsf + 1


C     Calculate chi2:
      chi2 = 0
      ndf = 0
      do if2=1,nf2tot
         f2 = f2vave(if2)
         ndf = ndf - 1
         if (mod(IAVMETH,10).eq.1) then
            fl = f2vave(if2+NF2TOT)
            if (fl.ne.0) then
               ndf = ndf - 1
            endif
         else
            fl = 0.        ! absorbed in "F2"
         endif
         
         do iexp=1,NMEAS
            if (F2TAB(if2,iexp).ne.0) then
               ndf = ndf + 1
               sum = F2TAB(if2,iexp)
               sred = f2 - fl * YFACTFL(if2,iexp)
               do isys=1,NSYSTOT
                  sum = sum + SYSTAB(isys,if2,iexp)*SYSSH(isys)
               enddo
               chi2loc = (sred-sum)
     $              /F2ETAB(if2,iexp)
               chi2 = chi2 + chi2loc**2
            endif
         enddo
      enddo
      print *,' '
      print *,' '
      print '(''chi2/dof (no systematics)   ='',F6.2,''/'',I6)'
     $     ,chi2,ndf
      do isys=1,NSYSTOT
         chi2 = chi2 + SYSSH(isys)**2
      enddo
      print '(''chi2/dof (with systematics) ='',F6.2,''/'',I6)'
     $     ,chi2,ndf
      
      print *,' '
      print *,' '

 17   format (F10.2,2F10.6,6F10.3)
 18   format (9A10)
      if (nsf.eq.2) then
         if (mod(IAVMETH,10).eq.1) then
            print 18,' q2 ',' x ',' y ',' F2 ', ' F2 stat', ' F2 tot'
     $           ,' FL ', ' FL stat', ' FL tot'
         else
            print 18,' q2 ',' x ',' y ',' F2 ', ' F2 stat', ' F2 tot'
     $           ,' xF3 ', ' xF3 stat', ' xF3 tot'
         endif
         print '(94(''-''))'
      else
         
      endif

      open(73,file='output/f2fl.dat',status='unknown')

      do if2=1,NF2TOT
         if (nsf.eq.2) then
            print 17,q2tab(if2),xtab(if2),ytabmax(if2)
     $           ,f2vave(if2)
     $           ,f2EstAve(if2)
     $           ,f2EAVE(if2)
     $           ,f2vave(if2+nf2tot)
     $           ,1/sqrt(diagst(if2+nf2tot))
C     $           ,f2EstAve(if2+nf2tot)
     $           ,f2EAVE(if2+nf2tot)
            write(73,17) q2tab(if2),xtab(if2),ytabmax(if2)
     $           ,f2vave(if2),f2EstAve(if2),f2EAVE(if2)
     $           ,f2vave(if2+nf2tot)
C     $           ,f2EstAve(if2+nf2tot)
     $           ,1/sqrt(diagst(if2+nf2tot))
     $           ,f2EAVE(if2+nf2tot)

C            print *,'ho',1/sqrt(diagst(if2+nf2tot)),
C     $           1/sqrt(diag(if2+nf2tot))

         else
         endif
      enddo

      close (73)
      print '(94(''-''))'

      print *,' '
      print *,' '
C---------------------------------------------------------
      end
