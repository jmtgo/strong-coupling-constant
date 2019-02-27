*
* Created 4 July 2004 by SG
* Include file for f2ave project
*
      integer nexpmax ! Max. number of experiments. Exp. share systematics
      integer nmeasmax! Max. number of measurements (H1 CC, NC etc)
      integer nf2max  ! max. number of F2 points
      integer nsystmax 
      integer NREACTMAX ! max. number of types of measurements
      parameter (NEXPMAX=50,NMEASMAX=50,
C>>> $     NF2MAX=5000,NSYSTMAX=150,NREACTMAX=10)
     $     NF2MAX=5000,NSYSTMAX=200,NREACTMAX=10)

C Input tables:
     
      character*80 EXPNAME(NEXPMAX) ! experiment name

      integer NREACT                ! number of type of meas.
      integer RLIST(NREACTMAX)      ! List of reactions
      integer RCOUNT(NREACTMAX)     ! Number of points for each reaction type
C 28 june 06:
      real    RPOLAR(NREACTMAX)      ! Polarization for a given reaction

      integer NSYSR(NREACTMAX)          ! Total number of syst. affecting reaction
      integer ISYSR(NSYSTMAX,NREACTMAX) ! List of systematics

      integer NMEAR(NREACTMAX)          ! Number of measurements for each reaction type
      integer IMEAR(NMEASMAX,NREACTMAX) ! Number of measurements for each reaction type

      integer IREACT(NMEASMAX)      ! reaction type
      real  YMAXF(NMEASMAX)
      real  POLARF(NMEASMAX)

      character*132 SYSNAM(NSYSTMAX)
      real COMEXP            !  Center of mass energy for average

C
C Average cross sections or separate SF:
C
C  IAVMETH = 0  -- average x-sections
C  IAVMETH = 1  -- separate F2 and FL
C  IAVMETH = 10 -- separate F2 and xF3  (not yet)
C  IAVMETH = 11 -- separate F2, FL and xF3
C
      integer IAVMETH

C---------------------------------------------------------------
C FL constraints:   (FL - Raverage/(Raverage+1) F2)^2/(Rerror/(Rerror+1))^2
C
      logical LFLLIMITS                ! impose or not extra constraint
      real RAVERAGE                    ! default 0.5
      real RERROR                      ! default 0.5 Assumes F2=1 ...

      Common/CFLCONST/RAVERAGE,RERROR,LFLLIMITS
C
C end of FL constraints
C--------------------------------------------------------------


      real SYSTAB(NSYSTMAX,NF2MAX,NMEASMAX) ! Table of all syst. dependences
      real SYSTABOrig(NSYSTMAX,NF2MAX,NMEASMAX) ! Table of all syst. dependences
      real F2TAB(NF2MAX,NEXPMAX)            ! Table of all meas. values 
      real F2ETAB(NF2MAX,NEXPMAX)           ! Table of all errors.
      real F2ETABOrig(NF2MAX,NEXPMAX)           ! Table of all original errors.

C Save unorrelated and stat. errors separately:
      real F2ETAB_STA(NF2MAX,NEXPMAX) 
      real F2ETAB_UNC(NF2MAX,NEXPMAX) 
      real F2ETAB_STAORIG(NF2MAX,NEXPMAX)

      integer REACTAB(NF2MAX)         ! reaction type for this point
      real XTAB(NF2MAX)               ! x  of a point
      real Q2TAB(NF2MAX)              ! q2 of a point
      real CMETAB(NF2MAX)             ! CME of a point

      real YTAB(NF2MAX)               ! y of a point -- taken from higer CME if correction is applied.
      real YTABORIG(NF2MAX,NEXPMAX)   ! keep also original y
      real YTABMAX(NF2MAX)            ! maximal Y
      real YFACTFL(NF2MAX,NEXPMAX)    ! for FL fit yfactors
      real YFACTF3(NF2MAX,NEXPMAX)    ! for xF3 fit yfactors


      integer NEXP,NMEAS,NF2(NMEASMAX),NSYST(NEXPMAX),NF2TOT,NSYSTOT,
     $     NSFSEC
      integer INSYST(2,NEXPMAX)  ! reference to first/last source

      integer NSYSM(NMEASMAX)    ! number of syst. for a measurement
      integer ISYSM(NSYSTMAX,NMEASMAX) ! reference from meas. syst --> syst
      real YMAXINT               ! Limit for x-sect. integration
      real YMAXAVE(NREACTMAX)    ! Limit for x-section measurement for different CMS
      
      integer NMEASF2(NF2MAX)    ! number of measurements per point
      
C Control:
      integer IDEBUG
      integer IOUT
      integer ICHECK


C Add extra stat. error adjustment:
      integer nstatadj

C R influence on average:
      integer IRSTUDY

C
C Treat an error source as absolute or as relative
C
      logical 
     $     relative(0:NSYSTMAX), userelative

C
C Sqrt scaling for stat. errors:
C       
      logical LSqrtStat,LSqrtStat2,LConstStat

C
C Log (S_red) average:
C
      logical LOGAVE

      logical LFractal

C To study biases, simulate XS:
      integer nsim,iseed
  
C Swim control:
      integer IPDG
      logical LHERA02

C Average control:
      logical LAveSameExp

      common/cdebug/idebug,iout
      common/CN/SYSNAM,EXPNAME
      Common/CF2/XTAB,Q2TAB,F2TAB,F2ETAB,F2ETABOrig,SYSTAB,COMEXP,
     $     SYSTABOrig,YTAB,CMETAB,YFACTFL,YFACTF3,YTABORIG,YTABMAX,
     $     F2ETAB_STA,F2ETAB_UNC,F2ETAB_STAORIG,NMEASF2
      Common/CF2I/NEXP,NMEAS,NF2,NSYST,NF2TOT,NSYSTOT,INSYST,IPDG,
     $     IREACT,NREACT,REACTAB,RLIST,RCOUNT,NSYSM,ISYSM,NSYSR,ISYSR,
     $     NMEAR,IMEAR,ICHECK,YMAXINT,RPOLAR,YMAXF,POLARF,YMAXAVE
     $     ,LFractal,IRSTUDY,IAVMETH,NSFSEC, LSqrtStat,LSqrtStat2
     $     ,LHERA02,LAveSameExp,LConstStat

C Absolute/relative errors:
      Common/CSTAT/nstatadj,relative,userelative,LOGAVE

C Simulation:
      Common/CSIM/NSIM,iseed

C Output:
      real F2Vave(NF2MAX)
     $     ,F2Eave(NF2Max)
     $     ,F2EstAve(NF2Max)
     $     ,F2EstAveTrue(NF2Max)
     $     ,F2Esyave(NF2Max)
      real SYSSH(NSYSTMAX)       ! Shifts of the systematics
      real*8 ERRSYST(NSYSTMAX)
      Common/CF2out/F2Vave,F2Eave,F2Estave,F2EstAveTrue,F2Esyave,syssh
      Common/CF2outd/ERRSYST

