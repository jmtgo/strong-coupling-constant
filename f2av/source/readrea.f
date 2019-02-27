      Subroutine readrea(IREA,N,NSYS,XQY,XS,XSE,SYS,CFILE)
C----------------------------------------------------------
C
C Open output/tab.dat file and read reaction data
C
C-----------------------------------------------------------
      implicit none
      
      character*(*) CFILE
      integer IREA,N,NSYS,irec
      real XS(N),XSE(N),XQY(3,N)
      real sys(100,N)

      open (51,file=CFILE,status='old',err=51)
      read (51,*)
      read (51,*)
 17   call readtable(n,nsys,irec,xqy,xs,xse,sys)
      if (n.gt.0 .and. irec.ne.irea) then
         n = 0
         goto 17
      endif
      if (n.gt.0) then
      else
      endif
      

      close (51)
      return
 51   print *,'Error opening tab.dat file. Stop'
      stop
      end
