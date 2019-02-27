      Subroutine readtable(N,NSYS,IREC,XQY,XS,XSE,SYS)
C
C Read a table in QCDFIT format
C
      implicit none

      integer N,IREC,NSYS
      real XS(N),XSE(N),comexp,XQY(3,N)
      real SYS(100,N)
      real tmp
      integer i,j,NSYStmp,isys,list(100)
      character*8 cnam,ctmp
C-----------------------------------------------------
 1717 format (i5,a8,2i3,i5,F10.3,1X,a1,100i2)
 1718 format (8E14.7)
 1719 format (200F8.4)

      read (51,1717,err=99,end=99) N,
     $     CNAM,
     $     NSYS,NSYStmp,IREC,
     $     COMEXP,ctmp,(list(isys),isys=1,NSYS)

      print *,NSYS,N,list(nsys)
C Read data points:
      do i=1,N
         read (51,1718) xqy(1,i),xqy(2,i),xqy(3,i),
     $       XS(i),tmp,XSE(i),tmp
C Translate to relative errors:
         XSE(i) = XSE(i)/XS(i)
         read (51,1719) (SYS(j,i),j=1,NSYS)
         do j=1,NSYS
C Translate to relative (from %)
            SYS(j,i) = SYS(j,i)/100.
         enddo

C         print *,xqy(1,i),xqy(3,i),xs(i),xse(i),tmp
      enddo

C-----------------------------------------------------
      return
      
 99   n = 0
      end
