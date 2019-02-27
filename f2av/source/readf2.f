      Subroutine readf2
C---------------------------------------------------------------
C
C Created 4 July 2004 by SG.
C Read F2 tables
C
C---------------------------------------------------------------      
      implicit none
      include '../include/f2.inc'
C Locals:
      integer NF2Files,NFilesMax
      parameter (NFilesMax=100)
      character*132 F2Files(NFilesMax),Q2XGRID


      real CMEYMAX
      real YRef
      Logical LSwimY  ! Swimming mimimizes change in Y, not in X
      Logical LCorrCmeAfter ! correct CME after the swimming 

      Logical LCORCMSE  ! 7 mar 05, correct CMS energy difference
      namelist/INPUT/NF2Files,NEXP,F2Files,Q2XGRID,IDEBUG,IOUT,IPDG
     $     ,EXPNAME,COMEXP,NREACT,RLIST,ICHECK,LCORCMSE,YMAXINT,RPOLAR
     $     ,NSTATADJ,NSIM,ISEED,YMAXAVE,LFRACTAL,IRSTUDY,IAVMETH
     $     ,LFLLIMITS,RAVERAGE,RERROR, LSqrtStat, LSqrtStat2,logave
     $     ,LHERA02,LAveSameExp,LConstStat,CMEYMAX, LSwimY
     $     ,LCorrCmeAfter

      namelist/SYST/NSYSTOT,NSYST,SYSNAM
     $     ,userelative,relative

      integer ifile
      character*2024 ctmp
      logical lfline
      integer nsystloc,nf2loc
      integer nqloc,nxloc
C Q2/X grid:
      integer nqmax,nxmax
      parameter (nqmax=1000,nxmax=1000)
      real q2buf(nqmax),xbuf(nxmax,nqmax)
      integer NCMEMAX
      parameter (NCMEMAX=5)

      integer nf2exp(NCMEMAX,NREACTMAX,nxmax,nqmax)
      integer if2exp(NCMEMAX,NREACTMAX,nxmax,nqmax)

      integer NCME,ICME
      real CMELIST(NCMEMAX)

      integer ninput

C File header:
      integer 
     $     nf2file,
     $     nsfile
      integer SYSTYPE(NSYSTMAX) ! ref. to systematics

      logical
     $     percent(0:NSYSTMAX)
      
      integer REACTION
      real CMEEXP             ! 7 Mar 05
      real YMAXEXP            ! 25 June 06
      real polar              ! 28 june 06
      real YSWIMEXP           ! VC. 18.11.2012
      namelist/FHEAD/nf2file,nsfile,systype,percent,
     $     reaction,cmeexp,
     $     ymaxexp,polar,YSWIMEXP

      integer ncorr
      real q2,x,f2,f2st,f2sys(nsystmax)   
      real Scale
      real y,CME,yorig

      real f2st_true, f2uncor

      integer IRea  ! for current data pointer to the reaction list

      integer iq2,ix,iexp
      real xd,q2d,xdmin,q2dmin
      integer i,j,if2,k,nnn
      real W1,W2

      real signF3

      integer systmp(NSYSTMAX),isub
C---------------------------------------------------------------      
      ninput = 0 
C Zero F2/F2E arrays:
      do i=1,NF2MAX
         do j=1,NEXPMAX
            F2TAB(i,j)  = 0.0
            F2ETAB(i,j) = 0.0
            F2ETAB_STA(i,j) = 0.0
            F2ETAB_STAORIG(i,j) = 0.0
            F2ETAB_UNC(i,j) = 0.0
            F2ETABOrig(i,j) = 0.0
         enddo
      enddo

C Defaults:
      IPDG = 5003
      NEXP = 0
      YMAXINT = 1.

      NSTATADJ = 0
      NSIM = 0
      ISEED = 17

      USERelative = .false. ! Ignore "relative" table

      do i=0,NSYSTMAX
         relative(i) = .false.
      enddo

      do i=1,NREACTMAX
         YMAXAVE(i) = 1.        ! average all data
      enddo

      LSqrtStat  = .false.       ! do not consider stat. errors separately.
      LSqrtStat2 = .false.       ! also do not rescale stat. errors using syst. shifts
      LConstStat = .false.
      LOGAVE     = .false.       ! average S_red, not log(S_red)
      LHERA02    = .false.

      LAveSameExp = .false.      ! do not average points from the same dataset

      LCorrCmeAfter = .false.

      LFractal = .false.
      IRSTUDY = 0
      IAVMETH = 0         ! by default simple average

      LFLLIMITS = .false. ! no extra term
      RAVERAGE  = 0.5
      RERROR    = 0.5

      NCME = 1

C 7 Mar 05:
      LCORCMSE = .false.  ! Do not correct CMS energy difference by def.

C 28 Feb 2010:
      CMEYMAX = 0.

C 28 Feb 2010:
      LSwimY = .false.

C Read steering cards:

      read (5,nml=INPUT,ERR=91,END=91)
      print '(''Read F2 from'',i4,'' files'')',NF2Files
      CMELIST(1) = COMEXP       ! At least one CME

C Select method for stat. errors correction:
      if (LSqrtStat.and.LConstStat) then
         print *,
     $    'ERRROR, INCONSISTENT STEERING, LSqrtStat.and.LConstStat=true'
         print *,'STOP'
         stop
      endif
      
         
C 10 Nov 2007: Do not correct for CME if SF F2 and FL are extracted
      if (mod(IAVMETH,10).eq.1) then  
         print *,' '
         print *,' Extract SF FL '
         print *,' Reset CME correction to none '
         LCORCMSE = .false.
         if (LFLLIMITS) then
            print *,' Impose soft FL limits'
            print *,' R mean=',RAVERAGE,' R sigma=',RERROR
         endif
         print *,' '
      endif

      if (LCORCMSE) then
         print *,'Correct CMS eneregy to ',COMEXP

         if (YMAXAVE(1).ne.1.) then
            print '(''Average only data with y<'',10F6.2)',
     $           (ymaxave(i),i=1,NREACT)
            print '(''for reactions            '',10I6)',
     $           (RLIST(i),i=1,NREACT)
         endif

         if (CMEYMAX.ne.0) then
            print '(''Use reference CME for Ymax='',F12.0)',CMEYMAX
         endif

      else
         print *,'Leave CMS energy per experiment unmodified'
      endif

      if (LFractal) then
         print *,'Use Fractal fit for CME/Swimming corrections'
      endif

      if (LSqrtStat) then
         print *,
     $ 'LSqrtStat = True,  rescale stat. and uncorr errors separately'
         if (LSqrtStat2) then
            print *,
     $     'LSqrtStat2 = true, correct syst. bias for stat. errors'
         else
            print *,
     $  'LSqrtStat2 = false, ignore syst. bias for stat. errors'
         endif
      else
         print *,'LSqrtStat = False'
      endif

      if (LConstStat) then
         print *,
     $ 'LConstStat = True,  keep Stat. errors fixed'
      else
         print *,'LConstStat = False'
      endif

      if (LHERA02) then
         print *,
     $      'Use HERAPDF0.2 parameterization for swimming corrections'
      endif

      if (LAveSameExp) then
         print *,
     $'Average  data from the same dataset swimming to the same x,Q^2'
      endif

C Read also systematic table namelist:
      read (5,nml=SYST,ERR=91,END=91) 

      print *,'Use relative syst. errors=',UseRelative

      if (UseRelative) then
         print '(A3,'' '',A30,A10)',' ','Source Name              '
     $        , 'Relative'
         do i=1,NSYSTOT
            print '(I3,'' '',A30,L10)', i,SYSNAM(i),Relative(i)
         enddo
      endif

      if (LSwimY) then
         print '(''Swim to closest Y point'')'
      else
         print '(''Swim to closest X point'')'
      endif

      if (NEXP.eq.0) then
C Each experiment is one input file:
         NEXP = NF2FILES
      endif

      NMEAS = NF2FILES ! one file --> one measurement

C Fill references experiment --> systematic
      do iexp=1,NEXP
         if (iexp.eq.1) then
            INSYST(1,iexp) = 1
         else
            INSYST(1,iexp) = INSYST(2,iexp-1) + 1
         endif
         INSYST(2,iexp) = INSYST(1,iexp) + NSYST(iexp) - 1
      enddo

C      print *,(INSYST(1,iexp),INSYST(2,iexp),iexp=1,nexp)

      print *,' '

C Read Q2/X grid:
      open (51,file=q2xgrid,status='old',err=99)
      lfline = .true.
 11   read (51,'(A)',err=98,end=12) ctmp         
      if (ctmp(1:1).ne.'*') then
         if (lfline) then
            read (ctmp,*) nqloc,nxloc
            print *,nqloc,nxloc
            lfline = .false.
            j = 1
         else
c            print *,ctmp
            read (ctmp,*) q2buf(j),(xbuf(i,j),i=1,nxloc)
            j = j + 1
         endif
      endif
      goto 11
 12   close (51)
      print *,'Q2 grid=',(q2buf(i),i=1,nqloc)

C Read individual files with measurements:
      do ifile = 1, NF2Files
         CTMP = F2FILES(IFILE)
         print *,' '
         print '(''Reading file '',A)',CTMP(1:index(ctmp,' '))
         open (51,File=F2Files(ifile),status='old',err=92)
C reset to defaults
         NSFILE = 0
         NF2FILE = 0
         do i=1,NSYSTMAX
            systype(i) = 0       ! assume un-correlated
            percent(i) = .false. ! assume absolute errors
            systmp(i) = 0
         enddo
         percent(0) = .false.    ! stat. error is not in percent
         CMEEXP     = COMEXP     ! By default the same CofM Energy
         YMAXEXP    = YMAXINT    ! By default the same Ymax
         POLAR      = 0.         ! No polarization by default
         YSWIMEXP   = 1.         ! No Q2-Y swimming for Y<YSWIMEXP, VC 18.11.2012
c read the header:
         read (51,nml=FHEAD,err=93,end=92)

         if (systype(1).eq.0) then
            isub = 1
         else
            isub = 0
         endif

         do i=1,nsfile
            if (systype(i).gt.0) then
               systmp(systype(i)) = i - isub
            endif
         enddo
         print '(''N meas, reaction, nsyst='',i5,2i5)'
     $        ,nf2file,reaction,nsfile
         print '(''Center of Mass energy='',F10.2)',CMEEXP


         
         print '(''SYSTFILE:'',A20,I5,50I4)',CTMP(1:index(ctmp,' ')),
     $        nsfile,(systype(i),i=1,nsfile)

         print '(''SYSFILET:'',A20,110I3)',CTMP(1:index(ctmp,' ')),
     $        (systmp(i),i=1,110)

         do i=1,NCME
            if (abs(CMELIST(i)-CMEEXP).lt.0.001) then
               goto 8712
            endif
         enddo
         print *,'Got new CME',CMEEXP
         NCME = NCME + 1
         CMELIST(NCME) = CMEEXP
 8712    continue


C for xF3 extraction, reset 515 to 615
         SignF3 = 1.
         if (IAVMETH.eq.10) then
            if (Reaction .eq. 515) then
               Reaction = 615
               SignF3   = -1.
            endif
         endif


C check if this reaction is on the reaction list:
         do j=1,NREACT
            if (Reaction.eq.RLIST(j)) then
               IREA = j
               goto 51
            endif
         enddo
C Error..
         print '(''Can not find reaction ='',i5)',Reaction
         print '(''Abort'')'
         Stop
 51      Continue

         IREACT(ifile) = REACTION

C Store also Ymax for integration:
         YMAXF(ifile) = YMAXEXP
C Store polarization:
         POLARF(ifile) = POLAR

C determine # of correlated systematic uncertainties:
         NCORR = 0
         do i=1,NSFILE
            if (systype(i).ne.0) then
                NCORR = NCORR + 1
                ISYSM(NCORR,ifile) = SYSTYPE(i)  ! Store reference
C 22 July 04: also store references for the given reaction:
                do j=1,NSYSR(irea)
                   if (ISYSR(j,irea).eq.SYSTYPE(i)) goto 32 ! Already present
                enddo
C Add new:
                NSYSR(IREA) = NSYSR(IREA) + 1
                ISYSR(NSYSR(IREA),IREA) = SYSTYPE(i)
 32             Continue
            endif
         enddo
         print '(''Correlated systematic='',i3)',ncorr
         NF2(ifile) = nf2file
         NSYSM(ifile) = NCORR

C 22 July 04: Reaction --> measurement reference
         NMEAR(IREA) = NMEAR(IREA) + 1
         IMEAR(NMEAR(IREA),IREA) = IFILE

         j = 0

 21      read (51,'(A)',err=93,end=22) ctmp         
         if (ctmp(1:1).ne.'*') then
c read the F2 points:
            j = j + 1
            if (j.gt.NF2file) then
               print 
     $              '(''Warning: too many F2 points'',i2,'' skip..'')',
     $                 j
               goto 22
            endif
            ninput = ninput + 1

            read (ctmp,*) q2,x,f2,f2st,(f2sys(i),i=1,nsfile)
C
C Convert %errors to absolute:
C
            if (percent(0)) then
               f2st = 0.01*f2st * f2
            endif
            do i=1,nsfile
               if (percent(i)) then
                  f2sys(i) = 0.01*f2sys(i) * f2
               endif
            enddo
C
C Be ready to average log of the x-section
C
            if (LOGAVE) then
               f2st = f2st/f2
               do i=1,nsfile
                  f2sys(i) = f2sys(i)/f2
               enddo
               f2 = log(f2)
            endif

C
C Store stat. and uncorrelated errors separately:
C
            f2st_true = f2st
            f2uncor   = 0. 
C
C Depending on systematic type add to statistics or not:
C
            do i=1,nsfile
               if (systype(i).eq.0) then
                  f2st     = sqrt(f2st**2+f2sys(i)**2)
                  f2uncor  = sqrt(f2uncor**2+f2sys(i)**2)
                  f2sys(i) = 0.
               endif
            enddo

c Check corresponding point on the grid:
            q2dmin = 1.E10
            do i=1,nqloc
C For now simple linear closest:
               q2d = abs(q2buf(i)-q2)
               if (q2d.lt.q2dmin) then
                  iq2 = i
                  q2dmin = q2d
               endif
            enddo



C 7 Mar 2005 SG: Correct for different CMS energy:
C 21 Apr 2007: Apply CME correction only if Y is smaller than YMAXAVE

            Y = Q2/(X*CMEEXP)
            Yorig = Y

C Now find X
C 28 Feb 2010, SG: Swim to closest y point, if generally requested            
C 18.11.2012,  VC: Swim to closest y point, if requested per file
            if (Q2.ne.q2buf(iq2).and.y.gt.YSWIMEXP) then
            print '(''Swim in Q2 at constant Y='',F6.4,
     $      '' to closest point in X='',F9.3,''/Y/'',F7.0,''='',F9.7)',
     $      y,q2buf(iq2),CMEEXP,q2buf(iq2)/(y*CMEEXP)
            endif

               xdmin = 1.E10
               do i=1,nxloc
C  Linear:
                if (LSwimY.or.(Q2.ne.q2buf(iq2).and.y.gt.YSWIMEXP)) then
                  xd = abs(q2buf(iq2)/(xbuf(i,iq2)*CMEEXP)-y) 
                else
                  xd = abs(xbuf(i,iq2)-x) 
                endif

                if (xd.lt.xdmin) then
                  ix = i
                  xdmin = xd
                endif
               enddo
C
C 28 Feb 2010: allow to calculate YRef using a fixed reference CME:
C
            if (CMEYMAX.eq.0) then
               YRef = Y
            else
               if (LCorrCmeAfter) then
                  YRef = Q2buf(iq2)/(Xbuf(ix,iq2)*CMEYMAX)
               else
                  YRef = Q2/(X*CMEYMAX)
               endif
            endif


            if (abs(CMEEXP-COMEXP).gt.0.001
     $           .and.LCORCMSE) then
               if (
     $              YRef.gt.YMAXAVE(IREA)
C  01.11.12 VC: do CME correction also for high y if CME difference is < 0.5%
     $              .and.abs(CMEEXP-COMEXP).gt.0.005*abs(COMEXP)
     $              ) then
                  CME = CMEEXP
                  print *,'Keep low CME y=',Yref,' YaveMax='
     $                 ,YMAXAVE(irea)
               else
                  Call XSCorrCME(Reaction,q2,x,CMEEXP,COMEXP,F2) 
                  CME = COMEXP
               endif
            else
               CME = CMEEXP                  
            endif

C Now find CME:
            ICME = 0
            do i=1,NCME
               if (abs(CMELIST(i)-CME).le.0.001) then
                  ICME = I
                  goto 1191
               endif
            enddo


            if (ICME.eq.0) then
               print *,'Problem finding CME, stop'
               stop
            endif


 1191       Continue

C Reset ICME to 1 for F2,FL fit:
            if (mod(IAVMETH,10).eq.1) then
               ICME = 1
            endif

C 5 June 2004 SG: We need to swim to this new central value ...
            Call F2Swim(Reaction,q2,x,q2buf(iq2),xbuf(ix,iq2),scale
     $           ,YMaxExp,Polar,CME)

            Y = Q2buf(iq2)/(Xbuf(ix,iq2)*CME)

            if (IDEBUG.ge.1) then
               print *,'Swim from, to:',q2,x,q2buf(iq2),xbuf(ix,iq2)
            endif


            if (LOGAVE) then
               F2 = F2 + log(scale)
            else
               F2 = F2*scale
               F2st = F2st*scale
               do i=1,nsfile
                  f2sys(i) = f2sys(i)*scale
               enddo
C Feb 5: Bug fix: scale original stat. error  too
               f2st_true = f2st_true * scale
               f2uncor   = f2uncor * scale
            endif

C Got X,Q2, check if it is a new one:
            nf2exp(icme,irea,ix,iq2) = nf2exp(icme,irea,ix,iq2) + 1
            
            if (nf2exp(icme,irea,ix,iq2).gt.1) then
               if2 = if2exp(icme,irea,ix,iq2)
            else
               if2 = NF2TOT + 1
               if2exp(icme,irea,ix,iq2) = if2               
               nf2TOT = if2
            endif
C Check if there is already a point at this X,Q2 for this experiment:
            if (F2TAB(if2,ifile).ne.0) then
               print *,'WARRNING: same X,Q2 twice for one experiment'
               print *,if2,xtab(if2),q2tab(if2),f2,f2tab(if2,ifile)

               if (LAveSameExp) then
                  
                  print *,'Orig:'
                  print *,'sigma:',F2TAB(if2,ifile),F2
                  print *,'stat :',F2ETAB_STA(if2,ifile),f2st_true
                  print *,'uncor:',F2ETAB_UNC(if2,ifile),f2uncor
                  print *,'totuc:',F2ETAB(if2,ifile),f2st

C SG: 8 June 2009 : Average data from the same experiment using stat. errors 
C
                  W1 = 1/F2ETAB_STA(if2,ifile)**2
                  W2 = 1/f2st_true**2
                  
                  W1 = W1/(1/F2ETAB_STA(if2,ifile)**2+1/f2st_true**2)
                  W2 = W2/(1/F2ETAB_STA(if2,ifile)**2+1/f2st_true**2)

                  F2TAB(if2,ifile)  = W1*F2TAB(if2,ifile) + W2*f2
                  F2ETAB_STA(if2,ifile) = sqrt(
     $                 1/(1/F2ETAB_STA(if2,ifile)**2+1/f2st_true**2))

                  F2ETAB_STAORIG(if2,ifile) =  F2ETAB_STA(if2,ifile)
                  F2ETAB_UNC(if2,ifile) = sqrt(
     $                 W1**2*F2ETAB_UNC(if2,ifile)**2+
     $                 W2**2*f2uncor**2)

                  F2ETAB(if2,ifile) = sqrt( F2ETAB_STA(if2,ifile)**2
     $                 + F2ETAB_UNC(if2,ifile)**2)
                  F2ETABOrig(if2,ifile) = F2ETAB(if2,ifile)

                  print *,'Weights:',W1,W2
                  print *,'Average:'
                  print *,'sigma:',F2TAB(if2,ifile)
                  print *,'stat :',F2ETAB_STA(if2,ifile)
                  print *,'uncor:',F2ETAB_UNC(if2,ifile)
                  print *,'totuc:',F2ETAB(if2,ifile)
                  print *,' '

                  do i=1,nsfile
                     if (SYSTYPE(i).ne.0) then
                        SYSTAB(SYSTYPE(i),if2,ifile) = 
     $                       W1*SYSTAB(SYSTYPE(i),if2,ifile)
     $                       +W2*f2sys(i)
                        SYSTABOrig(SYSTYPE(i),if2,ifile) = 
     $                       W1*SYSTABOrig(SYSTYPE(i),if2,ifile)
     $                       +W2*f2sys(i)
                     endif
                  enddo
                  goto 21
               endif

            endif

C Store ...
            XTAB(if2)  = xbuf(ix,iq2)
            Q2TAB(if2) = q2buf(iq2)     
            YTAB(if2) = Y


            CMETAB(if2) = CME
            REACTAB(if2) = REACTION
            F2TAB(if2,ifile) = f2
            F2ETAB(if2,ifile) = f2st
            F2ETABOrig(if2,ifile) = f2st

C store separate stat. and uncorr errors:

            F2ETAB_STA(if2,ifile)     = f2st_true
            F2ETAB_STAORIG(if2,ifile) = f2st_true
            F2ETAB_UNC(if2,ifile)     = f2uncor


C Store original Y:
            YTABORIG(if2,ifile) = Yorig 
            if (yorig.gt.YTABMAX(if2)) then
               YTABMAX(if2) = Yorig
            endif
C Add FL/xF3 factors:
            YFACTFL(if2,ifile) = Yorig**2/(1+(1-Yorig)**2)
            YFACTF3(if2,ifile) = (1-(1-Yorig)**2)
     $           /(1+(1-Yorig)**2)*SignF3
            


            do i=1,nsfile
               if (SYSTYPE(i).ne.0) then
                  SYSTAB(SYSTYPE(i),if2,ifile) = f2sys(i)
                  SYSTABOrig(SYSTYPE(i),if2,ifile) = f2sys(i)
               endif
            enddo
                        
         endif
         goto 21
         
 22      close (51)
      enddo

C Count how many points are of each reaction type
      do irea=1,NREACT
         do if2=1,NF2TOT
            if (REACTAB(if2).eq.RLIST(irea)) then
               RCOUNT(irea) = RCOUNT(irea) + 1
            endif
         enddo
      enddo

      print *,' '
      print '(''Read in '',i4,'' total  XS points'')',ninput
      print '(''Read in '',i4,'' unique XS points'')',nf2tot
      print *,' '

C Count in how many experiments each point is present
      do if2=1,NF2TOT
         NMEASF2(if2) = 0
         do iexp=1,NMEAS
            if (F2TAB(if2,iexp).ne.0) then
               NMEASF2(if2) = NMEASF2(if2) + 1
            endif
         enddo
      enddo

      if (idebug.ge.1) then
         print *,'All read F2:'
         do if2=1,NF2TOT
            print '(2F10.4,i5,50F8.2)'
     $           ,Q2TAB(if2),XTAB(if2),REACTAB(if2),
     $           (F2TAB(if2,i),i=1,NMEAS)
            do i=1,NMEAS
               if (F2TAB(if2,i).gt.0) then
                  nnn = nnn + 1
               endif
            enddo
         enddo
         print *,'Check total number of points=',nnn
      endif

      return
C---------------------------------------------------------------
 91   print *,'Error reading INPUT namelist'
      stop
 92   print *,'Can not open file:'
      print *,F2files(ifile)
      stop
 93   print *,'error reading file:'
      print *,F2files(ifile)
      stop
 98   print *,'error reading file:'
      print *,q2xgrid
      stop
 99   print *,'Can not open file:'
      print *,q2xgrid
      stop

      end
