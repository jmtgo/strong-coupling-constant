      Subroutine ReCovar(n,nsys,xsexp,xsth,xser,systab,var,covar,ifail)
C
C Re-calculate covariance matrix using different input 
C
      implicit none
      integer n,nsys
      real XSEXP(n),   ! experimental cross section
     $     XSTH(n),    ! theory input 
     $     XSER(n),    ! cross section errors (relative)
     $     SYSTAB(100,N)  ! cross section syst. errors (relative)

      
      real*8 VAR(1100,1100),COVAR(1100,1100) ! Variance/covariance matrix
      integer ifail

      real*8 TMP(1100)

      real*8 DIAG(1000),LAST(1100),CORR(100,1000),BOX(100,100)
      integer if2,isys,isys1,isys2,i,j
C------------------------------------------------------------------------
C Zero everything:
      do if2=1,N
         diag(if2) = 0.0
         do isys=1,Nsys
            corr(isys,if2) = 0.0
         enddo
      enddo

      do if2=1,N+NSYS
         last(if2) = 0.0
      enddo

      do i=1,N+NSYS
         do j=1,N+NSYS
            covar(i,j)  = 0
         enddo
      enddo
         

      do isys1=1,NSYS
         do isys2=1,NSYS
            box(isys1,isys2) = 0.0
         enddo
      enddo

      

C Diagonal:

      do if2=1,N
         diag(if2) = diag(if2) + 1/(XSER(if2)*XSEXP(if2))**2  ! use experiment here (?)
      enddo

C Last column:
      do if2=1,N
         last(if2) = 
     $              XSEXP(if2)/(XSEXP(if2)*XSER(if2))**2
      enddo

      do isys=1,NSYS
         do if2=1,N
            last(N+isys) = last(N+isys)
     $           - (SYSTAB(isys,if2)*XSTH(if2))
     $           * XSEXP(if2)/(XSEXP(if2)*XSER(if2))**2
         enddo
      enddo


C Correlation:
      do if2=1,N
         do isys=1,NSYS
            CORR(isys,if2) = CORR(isys,if2)
     $           - (SYSTAB(isys,if2)*XSTH(if2))
     $           /(XSEXP(if2)*XSER(if2))**2
         enddo
      enddo



C Box:

      do isys1=1,NSYS
         do isys2=1,NSYS
            do if2=1,N
               box(isys1,isys2) = box(isys1,isys2)
     $              + (SYSTAB(isys1,if2)*XSTH(if2))
     $              * (SYSTAB(isys2,if2)*XSTH(if2))
     $              / (XSEXP(if2)*XSER(if2))**2
            enddo
         enddo
      enddo

C Add 1 to box diagonal:

      do j=1,NSYS
         BOX(j,j) = BOX(j,j) + 1.0
      enddo


C Prepare CERN average:
C--------------------------------------------------------------

C Diagonal part:
      do if2=1,N
         covar(if2,if2) = diag(if2)
      enddo
C Correlation part:
      do if2=1,N
         do isys=1,NSYS
            covar(if2,N+isys) = corr(isys,if2)
            covar(N+isys,if2) = corr(isys,if2)
         enddo
      enddo
C "Box" part:
      do isys1=1,NSYS
         do isys2=1,NSYS
            covar(N+isys1,N+isys2) = box(isys1,isys2)
         enddo
      enddo
      
      Call DINV(n+nsys,covar,1100,tmp,ifail)

      if (ifail.ne.0) then
         print *,'Failed with covar ifail=',ifail
         stop
      endif

      do i=1,N+NSYS
         do j=1,N+NSYS
            var(i,j) = covar(i,j)  
         enddo
      enddo

C SG: invert reduced matrix
      Call DINV(n,var,1100,tmp,ifail)

      if (ifail.ne.0) then
         print *,'Failed with var ifail=',ifail
         stop
      endif

C------------------------------------------------------------------------
      end
