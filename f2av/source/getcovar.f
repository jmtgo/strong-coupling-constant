      Subroutine GetCovar(diag,last,corr,box,dfl)
C------------------------------------------------------------------------
C  Created 7 July 2004 by SG.
C  Calculate variance and co-variance matrix by brute force (dinv)
C
C------------------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
C Intermidiate dimensions:
      real*8 DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8 dfl(NF2MAX)
      real*8 Box(NSYSTMAX,NSYSTMAX)
C local:
      real*8 amat(nf2max+nsystmax,nf2max+nsystmax),tmp(nf2max+nsystmax)

      integer if2,isys,isys1,isys2,ifail,if2a,if2b
      real*8 errf2,errfl,corrf2fl
C------------------------------------------------------------------------

C Fill AMAT:

C Diagonal part:
      do if2=1,NSFSEC
         amat(if2,if2) = diag(if2)
      enddo

C FL calculation:
      if (mod(IAVMETH,10).eq.1) then
         do if2=1,NF2TOT
            amat(if2,if2+NF2TOT) = dfl(if2)
            amat(if2+NF2TOT,if2) = dfl(if2)
         enddo
      endif

C Correlation part:
      do if2=1,NSFSEC
         do isys=1,NSYSTOT
            amat(if2,NSFSEC +isys) = corr(isys,if2)
            amat(NSFSEC+isys,if2) = corr(isys,if2)
         enddo
      enddo
C "Box" part:
      do isys1=1,NSYSTOT
         do isys2=1,NSYSTOT
            amat(NSFSEC+isys1,NSFSEC+isys2) = box(isys1,isys2)
         enddo
      enddo
C Invert ...
      print *,'Get covariance matrix'
      print *,'Inverting...'
      Call DINV(NSFSEC+NSYSTOT,amat,nf2max+nsystmax,tmp,ifail)
      print *,'Done IFail=',ifail

C Write out covariance matrix:
      open (51,file='output/covar.dat',status='unknown')
C First dump x,q2:
      write (51,'(i6)') NSFSEC
      do if2=1,NSFSEC
         write (51,'(2i5,2F12.4)') if2,reactab(if2),q2tab(if2),xtab(if2)
      enddo
C Now dump correlation coefficients
      do if2a=1,NSFSEC
         write (51,'(i5,5000F6.2)') if2a,
     $        (amat(if2a,if2b)/sqrt(amat(if2a,if2a)*amat(if2b,if2b))
     $        ,if2b=1,NSFSEC)
      enddo
      close (51)

      if (IAVMETH.eq.1 .or. IAVMETH.eq.10) then
         open (51,file='output/f2flcor.dat',status='unknown')
         do if2=1,NF2tot
            errf2 = sqrt(amat(if2,if2))
            errfl = sqrt(amat(if2+NF2TOT,if2+NF2TOT))
            corrf2fl = amat(if2,if2+NF2TOT)/(errf2*errfl)
            write (51,'(F10.2,E12.4,3F10.3)') 
     $           q2tab(if2),xtab(if2)
     $           ,errf2
     $           ,errfl
     $           ,corrf2fl
         enddo
         close (51)
      endif
C
C Also, write out variance matrix !
C
      print *,'Get variance matrix'
      print *,'Inverting...'

C SG: invert reduced matrix
      Call DINV(NSFSEC,amat,nf2max+nsystmax,tmp,ifail)
      print *,'Done IFail=',ifail
C Write out covariance matrix:
      open (51,file='output/var.dat',status='unknown')
      write (51,'(i6)') NSFSEC
C First dump x,q2:
      if (mod(IAVMETH,10).eq.1) then
         do if2=1,NF2MAX
            write (51,'(2i8,3F18.10)') if2,reactab(if2),
     $           q2tab(if2),xtab(if2),f2vave(if2),f2vave(if2+NF2MAX)
         enddo

      else
         do if2=1,NSFSEC
            write (51,'(2i8,3F18.10)') if2,reactab(if2),
     $           q2tab(if2),xtab(if2),f2vave(if2)
         enddo
      endif
C Now dump correlation coefficients
      do if2a=1,NSFSEC
         write (51,'(i5,5000E18.10)') if2a,
     $        (amat(if2a,if2b),if2b=1,NSFSEC)
      enddo
      close (51)

C Diagonal correlation
      if (mod(IAVMETH,10).eq.1) then
         do if2a=1,NSFSEC
            do if2b=1,NSFSEC
               amat(if2a,if2b) = 0.
            enddo
         enddo

C     Diagonal part:
         do if2=1,NSFSEC
            amat(if2,if2) = diag(if2)
         enddo
         
C FL calculation:
         do if2=1,NF2TOT
            amat(if2,if2+NF2TOT) = dfl(if2)
            amat(if2+NF2TOT,if2) = dfl(if2)
         enddo
         
         Call DINV(NSFSEC,amat,nf2max+nsystmax,tmp,ifail)

         open (51,file='output/f2flcorst.dat',status='unknown')
         do if2=1,NF2tot
            errf2 = sqrt(amat(if2,if2))
            errfl = sqrt(amat(if2+NF2TOT,if2+NF2TOT))
            corrf2fl = amat(if2,if2+NF2TOT)/(errf2*errfl)
            write (51,'(F10.2,E12.4,3F10.3)') 
     $           q2tab(if2),xtab(if2)
     $           ,errf2
     $           ,errfl
     $           ,corrf2fl
         enddo
         close (51)
      endif

      end
