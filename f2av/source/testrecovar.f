      Program testrec
C------------------------------------------------------------------
C
C------------------------------------------------------------------
      implicit none

      integer nsys,n,i,ifail
      real xsexp(1000),xsth(1000),xser(1000),XQY(3,1000)
      real systab(100,1000)
      real*8 var(1100,1100),covar(1100,1100)

C-----------------------------------------------------------------
      n = 1
C read in data for IREA=615
      call readrea(615,N,NSYS,XQY,XSEXP,XSER,SYSTAB,'output/tab.dat')
      print *,n

      do i=1,1
         call recovar(n,nsys,xsexp,xsexp,xser,systab,var,covar,ifail)
      enddo
      print *,var(1,1),var(2,1),var(1,2),var(n,n)
      print *,'done'
      end
