      Subroutine Chi2Res(IOpt,N,XSth,X,Y,Q2,IReact,XSexp,Chi2,CFILE
     $     )
C
C Created 21 Apr 2007 by SG.
C Read the syst. source file CFILE and calulate Chi2
C for a given theoretical input XSth.
C
C Iopt = 0 --- Initialization: read CFILE, return
C           N -- N cross section points
C           X(N),Q2(N),XSexp(N) -- x,q2,type, and exp. value for
C           each measurement. XSth and Chi2 are undefined.
C Iopt = 1  Report N,X(N),Q2(N),,XSexp(N) without reading the file
C IOpt = 2  Calculate Chi2 using XSth as input
C
      implicit none
C I/O:
      integer IOpt,N
      real*8 XSth(N),X(N),Y(N),Q2(N),XSexp(N),Chi2
      
      character*(*) CFILE
      integer IReact

C Internal:
      
      integer NPMAX,NFMAX,ifail,NSYSMAX
      parameter (NPMAX=1000)
      parameter (NFMAX=10)
      parameter (NSYSMAX=100)

      character*256 cfiles(nfmax)

      integer NN(NFMAX),IR(NPMAX,NFMAX),nfiles,NSYS(NFMAX)

      real XS(NPMAX,NFMAX),XSE(NPMAX,NFMAX),SYS(100,NPMAX,NFMAX)
      real XYQ(3,NPMAX,NFMAX),XSTHs(NPMAX)

      real*8 VAR(NPMAX+NSYSMAX,NPMAX+NSYSMAX)
      real*8 COVAR(NPMAX+NSYSMAX,NPMAX+NSYSMAX)
      
      integer i,j,idum,ifile
      
      logical LInit
      data LInit/.false./

      data nfiles/0/
      save NN,IR,XS,VAR,LInit,nfiles,XYQ,NSYS
C-------------------------------------------------

      Chi2 = 0.

C Loop over variance matrix files

      do i=1,nfiles
         if (CFILE.eq.CFILES(i)) then
            ifile = i
            goto 89
         endif
      enddo
C Found new file:
      if (iopt.eq.0) then
         nfiles = nfiles + 1
         ifile = nfiles
         CFILES(nfiles) = CFILE
      else
C Error, new file but no initialization:
         print *,'No initialization call for variance matrix'
         print *,cfile
         print *,'Stop'
         stop
      endif

 89   continue

      if (IOpt.eq.0) then

         call readrea(ireact,NN(ifile),NSYS(ifile),XYQ(1,1,ifile),
     $        XS(1,ifile),XSE(1,ifile),SYS(1,1,ifile),CFILE)


         PRINT *,nsys(ifile),'bla'
C--- Initialization ---
         LInit = .true.
      endif

      if (.not. LInit) then
         print *,'Error in CHI2VAR --- called before initialization'
         print *,'STOP'
         stop
      endif

      if (IOpt.eq.0 .or. IOpt.eq.1) then
C--- Report the grid ---
         N = NN(ifile)
         do i=1,NN(ifile)
            X(i)     = XYQ(1,i,ifile)
            Y(i)     = XYQ(2,i,ifile)
            Q2(i)    = XYQ(3,i,ifile)
            XSExp(i) = XS(i,ifile)
         enddo
      elseif (IOpt.eq.2) then

C Re-calculate variance matrix:
         do i=1,NN(ifile)
            xsths(i) = Sngl(xsth(i))
         enddo
         Call ReCovar(NN(ifile),NSYS(ifile),xs(1,ifile),xsths, ! xs(1,ifile),
     $        xse(1,ifile),sys(1,1,ifile),var,covar,ifail)


C--- Calculate Chi2 ---
         do i=1,NN(ifile)
            do j=1,NN(ifile)
               Chi2 = Chi2 + 
     $              (XSth(i)-Dble(XS(i,ifile)))
     $              *(XSth(j)-Dble(XS(j,ifile)))*VAR(i,j)
            enddo
         enddo
      endif


      Return
 91   continue
      print *,'Failed to open Chi2 file',cfile
      stop
 92   continue
      print *,'Failed to read Chi2 file',cfile
      stop
      end
