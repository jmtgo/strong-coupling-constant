      Subroutine StatRecalc
C--------------------------------------------
C 
C Created 28 Mar 2007 by SG.
C
C Rescale stat. error of the measurement assuming that relative 
C errors do not change.
C--------------------------------------------
      implicit none
      include '../include/f2.inc'
      integer if2,iexp,isys
      real f2orig, f2new, erorig, scal
      real uncor, stat
C 9 jan 09: add also expected value - systematic shift:
      real f2newSys
C--------------------------------------------
      do if2=1,NF2TOT
         do iexp=1,NMEAS

            f2orig = F2TAB(if2,iexp)
            erorig = F2ETABOrig(if2,iexp)

            if (erorig.eq.0) goto 1717
C
C Special case: log ave and stat. errors
C
            if ( LOGAVE ) then

               if (.true.
     $              .and.NSTATADJ.lt.0
     $              .and.LSQRTSTAT) then
                  f2orig = exp(f2orig)
                  f2new  = exp(F2VAVE(if2))


                  stat  = F2ETAB_STAORIG(if2,iexp) * f2orig
                  scal  = sqrt( f2new/f2orig)
                  stat  = stat * scal/f2new 

                  uncor = F2ETAB_UNC(if2,iexp)
                  F2ETAB_STA(if2,iexp) = stat
                  F2ETAB(if2,iexp) = sqrt(stat**2+uncor**2)
               
C                  print *,'hi',f2orig,f2new,stat,uncor
               endif
C     Skip the rest (to the next point):
               goto 1717
            endif

C----
            if (NSTATADJ.gt.0) then
               f2new = f2orig
               do isys =1,NSYSTOT 
                  f2new = f2new + 
     $                 SYSTAB(isys,if2,iexp)*SYSSH(isys)
               enddo
            else
               f2new = F2VAVE(if2)
C 9 jan 09:
               f2newSys = f2new
               do isys=1,NSYSTOT
                  f2newSys = f2newSys -
     $                 SYSTAB(isys,if2,iexp)*SYSSH(isys)
               enddo
            endif

            if (f2orig.eq.0) then
               scal = 0.
            else
               scal = f2new/f2orig
            endif

C New error:
            if (LSqrtStat .or. LConstStat) then ! sqrt or no scaling

C Sqrt scaling for stat:
               if (LSqrtStat) then
                  if (LSqrtStat2) then
                     stat  = F2ETAB_STAORIG(if2,iexp) * 
     $                    sqrt(f2newSys/f2orig)
                  else
                     stat  = F2ETAB_STAORIG(if2,iexp) * sqrt(scal)
                  endif
               else
                  stat =  F2ETAB_STAORIG(if2,iexp)
               endif
               uncor = F2ETAB_UNC(if2,iexp)
               F2ETAB_STA(if2,iexp) = stat

C               print *,'old, new stat=',stat,F2ETAB_STAORIG(if2,iexp)

C Rescale uncor syst. part:
               if (.not.USERelative .or. Relative(0)) then
                  uncor = uncor * scal            
               endif

C Total uncor error:
               F2ETAB(if2,iexp) = sqrt(stat**2+uncor**2)
            else
               if (.not.USERelative .or. Relative(0)) then
                  F2ETAB(if2,iexp) = erorig * scal            
               endif
            endif

C New systematic sensitivity:
            do isys=1,NSYSTOT
               if (.not.USERelative .or. Relative(isys)) then
                  SYSTAB(isys,if2,iexp) = SYSTABOrig(isys,if2,iexp)
     $              * scal
               endif
            enddo
 1717    enddo
      enddo
C      print *,'ha',f2etab(1,1),f2etab(1,2),f2etaborig(1,1)
C     $     ,f2etaborig(1,2)
C      print *,'ho',systab(6,1,1),systab(16,1,2),systaborig(6,1,1)
C     $     ,systaborig(16,1,2)
C      print *,(systab(isys,1,2),isys=1,nsystot)
      end
