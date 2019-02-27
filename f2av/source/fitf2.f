      Subroutine fitf2
C---------------------------------------------------------------------------
C
C Created 4 July 2004 by SG.
C   Average F2 from different experiments
C
C---------------------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
C Intermediate dimensions:
      real*8 DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8 Box(NSYSTMAX,NSYSTMAX)
      real*8 DFL(NF2MAX), DXF3(NF2MAX),dflst(NF2MAX)
      real*8 DiagST(NF2MAX)

C A copy:
      real*8 boxs(NSYSTMAX,NSYSTMAX)
      real*8 Diags(NF2MAX),Lasts(NF2MAX),CORRs(NSYSTMAX,NF2MAX)

C--------------------------------------------------------------
 17   continue
      Call FillArrays(diag,last,corr,box,dfl,dxf3,diagst,dflst)      ! Prepare the system of equations

C Save a copy:
      call ucopy(box,boxs,nsystmax*nsystmax*2)
      call ucopy(diag,diags,nf2max*2)
      call ucopy(last,lasts,nf2max*2)
      call ucopy(corr,corrs,nsystmax*nf2max*2)

C Find F2 and systematics: 
      Call ToBlockDiag(diag,last,corr,box,dfl,dxf3,diagst,dflst)     ! Perform fast inversion


C If needed, recalculate stat errors and repeat the average:
      if (NSTATADJ.ne.0) then
         Call StatRecalc
         if (NSTATADJ.gt.0) then
            NSTATADJ = NSTATADJ - 1
         else
            NSTATADJ = NSTATADJ + 1
         endif
         goto 17
      endif

      if (mod(IOUT/10,10).eq.1) then
C We want to calculate covariance matrix 
         Call GetCovar(diags,lasts,corrs,boxs,dfl)     ! If requested, calculate covariance matrix
      endif

C--------------------------------------------------------------
      end
