      Subroutine Chi2Var(IOpt,N,XSth,X,Y,Q2,IRea,XSexp,Chi2,CFILE
     $     )
C
C Created 11 June 2006 by SG.
C Read the variance matrix stored in file CFILE and calulate Chi2
C for a given theoretical input XSth.
C
C Iopt = 0 --- Initialization: read CFILE, return
C           N -- N cross section points
C           X(N),Q2(N),IRea(N),XSexp(N) -- x,q2,type, and exp. value for
C           each measurement. XSth and Chi2 are undefined.
C Iopt = 1  Report N,X(N),Q2(N),IRea(N),XSexp(N) without reading the file
C IOpt = 2  Calculate Chi2 using XSth as input
C
      implicit none
C I/O:
      integer IOpt,N
      real*8 XSth(N),X(N),Y(N),Q2(N),XSexp(N),Chi2
      character*(*) CFILE
      integer IRea(N)

C Internal:
      
      integer NPMAX,NFMAX
      parameter (NPMAX=5000)
      parameter (NFMAX=10)

      character*256 cfiles(nfmax)

      integer NN(NFMAX),IR(NPMAX,NFMAX),nfiles
      real*8 XX(NPMAX,NFMAX),QQ(NPMAX,NFMAX),
     $     YY(NPMAX,NFMAX),XS(NPMAX,NFMAX)
      real*8 VAR(NPMAX,NPMAX,NFMAX)
      
      integer i,j,idum,ifile
      
      logical LInit
      data LInit/.false./

      data nfiles/0/
      save NN,IR,XX,YY,QQ,XS,VAR,LInit,nfiles
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
C--- Initialization ---
         LInit = .true.
         open (51,file=CFILE,status='OLD',err=91)
C Read in:
         read (51,'(i6)',err=92) NN(ifile)
         do i=1,NN(ifile)
            read (51,'(2i8,3F18.10)',err=92) idum,IR(i,ifile),
     $           qq(i,ifile),xx(i,ifile),xs(i,ifile)
         enddo
         do i=1,NN(ifile)
            read (51,'(i5,5000E18.10)',err=92) idum,
     $           (var(i,j,ifile),j=1,NN(ifile))
         enddo
         close (51)
         print '(''CHI2VAR: READ'',I4,'' CROSS SECTION POINTS'')',
     $        NN(ifile)
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
            X(i)     = XX(i,ifile)
            Q2(i)    = QQ(i,ifile)
            IRea(i)  = IR(i,ifile)
            XSExp(i) = XS(i,ifile)
         enddo
      elseif (IOpt.eq.2) then
C--- Calculate Chi2 ---
         do i=1,NN(ifile)
            do j=1,NN(ifile)
               Chi2 = Chi2 + 
     $              (XSth(i)-XS(i,ifile))
     $              *(XSth(j)-XS(j,ifile))*VAR(i,j,ifile)
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
