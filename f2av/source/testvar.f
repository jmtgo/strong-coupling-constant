      Program TestVar
C
C Created 11 June 2006 by SG. Test chi2var
C
      implicit none
C
      integer N,i
      real*8 XSTh(5000),X(5000),Y(5000),Q2(5000),XSEXP(5000),Chi2
      integer IREA(5000)
C   --------------------------------------------
      call Chi2Var(0,n,XSth,X,Y,Q2,IRea,XSexp,Chi2,'output/var.dat')
      do i=1,N
         print *,x(i),q2(i),irea(i),xsexp(i)
      enddo
C get a chi2:
      call Chi2Var(2,n,XSexp,X,Y,Q2,IRea,XSexp,Chi2,'output/var.dat')
      print *,chi2
      do i=1,1000
         XSexp(1) = XSexp(1)*0.99
         call Chi2Var(2,n,XSexp,X,Y,Q2,IRea,XSexp,Chi2,'output/var.dat')
         print *,i,xsexp(1),chi2
      enddo
      end
