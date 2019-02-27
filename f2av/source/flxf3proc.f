      Subroutine FLXF3PROC(diag,last,corr,corrT,dfl,dxf3,diagst,dflst)
C--------------------------------------------------------------------
C
C Created 10 Nov 2007 by SG.
C
C Pre-processing of the matrices for FL/xF3 case:
C-------------------------------------------------------------------
C  \   \            \
C   \   \            \   F2
C    \   \            \
C \   \       --->     \
C  \   \                \  FL
C   \   \                \
C-------------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 dfl(NF2MAX),dxf3(NF2MAX),CorrT(NSYSTMAX,NF2MAX),
     $     dflst(NF2MAX)
      real*8 diagst(NF2MAX)

      real*8 diagS(NF2MAX),lastS(NF2MAX)
      real*8 diagSS(NF2MAX)

      integer if2,isys,noff
      real*8 coef,coefst
C---------------------------------------------------------------------
      
      noff = NF2TOT

* Copy diag, last
      do if2=1,NSFSEC
         diagS(if2) = diag(if2)
         diagSS(if2) = diagst(if2)
         lastS(if2) = last(if2)
c         print *,'iaia',diag(if2),diagst(if2)
      enddo

      if (IAVMETH.eq.1  .or. IAVMETH.eq.10) then
* FL / xF3 separate pass:

* First diagonalise FL / xF3 part:
         do if2=1,NF2TOT
            coef = - DFL(if2)/DIAG(if2)
            DIAG(if2+noff) = DIAG(if2+noff) + coef * DFL(if2)

            coefst = - DFLST(if2)/DIAGST(if2)
            DIAGst(if2+noff) = DIAGst(if2+noff) + coefst * DFLST(if2)

c            print *,'hihi',DIAG(if2+noff), DIAGst(if2+noff)


            if (abs(DIAG(if2+noff)).lt.1.0D-10
     $           .or. YTABMAX(if2).lt.YMAXAVE(REACTAB(if2))) then
C This point has no sensitivity to FL (too low Y or single CME). Reset ...
c               print *,'reset'

               DIAG(if2+noff) = 0.0001
               DIAGst(if2+noff) = 0.0001
               Last(if2+noff) = 0.
               diagS(if2+noff) = 1.
               diagSS(if2+noff) = 1.
               lastS(if2+noff) = 0.
               DFL(if2) = 0.
               DFLST(if2) = 0.
               do isys=1,NSYSTOT
                  Corr(isys,if2+noff)  = 0.
                  CorrT(isys,if2+noff) = 0.
               enddo
            else

               Last(if2+noff) = Last(if2+noff) + coef * Last(if2)
               do isys=1,NSYSTOT
                  CorrT(isys,if2+noff) = CorrT(isys,if2+noff) 
     $                 + coef * CorrT(isys,if2)
               enddo
            endif
         enddo
C Now F2 part:
         do if2=1,NF2TOT
            coef = - DFL(if2)/diagS(if2+noff)
            DIAG(if2) = DIAG(if2) + coef * DFL(if2)


            coefst = - DFLST(if2)/diagSS(if2+noff)
            DIAGst(if2) = DIAGst(if2) + coef * DFLST(if2)

            Last(if2) = Last(if2) + coef * LastS(if2+noff)
            do isys=1,NSYSTOT
               CorrT(isys,if2) = CorrT(isys,if2) 
     $              + coef * Corr(isys,if2+noff)
            enddo

         enddo

         noff = noff + NF2TOT
      endif


      if (IDEBUG.gt.3) then
         print *,'syst. chec:'
         do if2=1,NSFSEC
            do isys=1,NSYSTOT
               print *,if2,isys,corrT(isys,if2),corr(isys,if2)
            enddo
         enddo
      endif
C No xF3 so far...

      end

