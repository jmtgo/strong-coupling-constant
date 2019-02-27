      Subroutine OutF2
C------------------------------------------------------------
C
C Created 7 June 04 by SG.
C
C Write out some useful output files.
C
C------------------------------------------------------------
      implicit none
      include '../include/f2.inc'
C Locals:
      real F2ETOT(NF2MAX,NEXPMAX)    ! Original total uncertainties
     $     ,F2ETOTN(NF2MAX,NEXPMAX)  ! New total uncertainties
      integer if2,if22,i,isys,iexp,nq2loc,j,irea
      logical F2FLAG(NF2MAX)
      real q2loc
      real sum(NEXPMAX)

      real F2tmp(NF2MAX),F2tmpe(NF2MAX)

      character*80 ctmp,crea
C------------------------------------------------------------


      if (LOGAVE) then
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               F2TAB(if2,iexp) = exp(F2tab(if2,iexp))
               F2ETAB(if2,iexp) = F2ETAB(if2,iexp)*F2TAB(if2,iexp)
               do isys=1,NSYSTOT
                  SYSTAB(isys,if2,iexp) = SYSTAB(isys,if2,iexp)
     $                 *F2TAB(if2,iexp)
               enddo
            enddo
         enddo
      endif

      if (mod(IOUT,10).eq.1       ! unmodified F2
     $    .or. mod(IOUT,10).eq.2  ! F2 shifted by systematics
     $    .or. mod(IOUT,10).eq.3  ! unmodified F2 y
     $    .or. mod(IOUT,10).eq.4  ! F2 shifted by systematics y
     $     ) then
C Create X,F2 files for each Q2 bin, write them into "output" directory
         do if2=1,NF2TOT
            F2FlAG(if2) = .false.
         enddo

C Get input total errors:
         do if2=1,NF2TOT
            do iexp=1,NMEAS
               F2ETOT(if2,iexp)  = F2ETAB(if2,iexp)**2
               F2ETOTN(if2,iexp) = F2ETAB(if2,iexp)**2
               do isys=1,NSYSTOT
                  F2ETOT(if2,iexp) = F2ETOT(if2,iexp) 
     $                 + SYSTAB(isys,if2,iexp)**2
                  F2ETOTN(if2,iexp) = F2ETOTN(if2,iexp) 
     $                 + SYSTAB(isys,if2,iexp)**2*ERRSYST(isys)**2
               enddo
               F2ETOT(if2,iexp) = sqrt(F2ETOT(if2,iexp))
               F2ETOTN(if2,iexp) = sqrt(F2ETOTN(if2,iexp))
           enddo
         enddo 

C Over reactions:
         do irea=1,NREACT
            write (crea,'(I10)') RLIST(irea)            
            i = 1
            do while ( crea(i:i) .eq. ' ')
               i = i + 1
            enddo
            crea = crea(i:80)

C Q2 info:
            ctmp = 'output/'//crea(1:index(crea,' ')-1)//'q2.dat'
            open (52, file=ctmp, status='unknown')

C
C Add a table with all Q2 data:
C
            ctmp = 'output/'//crea(1:index(crea,' ')-1)//'.dat'
            open (50,file=ctmp,status='unknown')
            
C Scan over all F2 points:         
            do if2=1,NF2TOT
               if (REACTAB(if2).ne.RLIST(irea)) goto 99

               if (.not.F2FLAG(if2)) then
                  F2FLAG(if2) = .true.
                  Q2loc = Q2TAB(if2)
                  nq2loc = 1
C Open the new Q2 file:
                  write (ctmp,'(F8.2,''.dat'')') q2loc
                  i = 1
                  do while ( ctmp(i:i) .eq. ' ')
                     i = i + 1
                  enddo
                  ctmp = ctmp(i:80)
                  ctmp = 'output/'//
     $                 crea(1:index(crea,' ')-1)//'_'
     $                 //ctmp(1:index(ctmp,' '))
C     print *,ctmp
                  open (51,file=ctmp,status='unknown')
 17               format (F13.8,20E13.5)
 117              format (2E13.6,20E13.5)
                  if (mod(IOUT,10).eq.1) then
                     j = 0
                     do iexp=1,NMEAS
                        if (RLIST(irea).eq.IREACT(iexp)) then
                           j = j + 1
                           f2tmp(j)  = f2tab(if2,iexp)
                           f2tmpe(j) = f2etot(if2,iexp)
                        endif
                     enddo
                     write (51,17) xtab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)
C All Q2 info:
                     write (50,117) q2loc,xtab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)
                  elseif (mod(IOUT,10).eq.2) then
C Float each point according to syst:
                     j = 0
                     do iexp=1,NMEAS
                        if (RLIST(irea).eq.IREACT(iexp)) then
                           j = j + 1
                           f2tmp(j)  = f2tab(if2,iexp)
                           f2tmpe(j) = f2etotn(if2,iexp)
                           do isys=1,NSYSTOT
                              f2tmp(j) = f2tmp(j) + 
     $                             SYSTAB(isys,if2,iexp)*SYSSH(isys)
                           enddo

                        endif
                     enddo

                     write (51,17) xtab(if2)
     $                    ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)

C All Q2 info:
                     write (50,117) q2loc,xtab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)



C add Y:
                  elseif (mod(IOUT,10).eq.3) then
                     j = 0
                     do iexp=1,NMEAS
                        if (RLIST(irea).eq.IREACT(iexp)) then
                           j = j + 1
                           f2tmp(j)  = f2tab(if2,iexp)
                           f2tmpe(j) = f2etot(if2,iexp)
                        endif
                     enddo
                     write (51,17) xtab(if2),ytab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)
C All Q2 info:
                     write (50,117) q2loc,xtab(if2),ytab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)

C add Y:
                  elseif (mod(IOUT,10).eq.4) then
C Float each point according to syst:
                     j = 0
                     do iexp=1,NMEAS
                        if (RLIST(irea).eq.IREACT(iexp)) then
                           j = j + 1
                           f2tmp(j)  = f2tab(if2,iexp)
                           f2tmpe(j) = f2etotn(if2,iexp)
                           do isys=1,NSYSTOT
                              f2tmp(j) = f2tmp(j) + 
     $                             SYSTAB(isys,if2,iexp)*SYSSH(isys)
                           enddo

                        endif
                     enddo

                     write (51,17) xtab(if2),ytab(if2)
     $                    ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)

C All Q2 info:
                     write (50,117) q2loc,xtab(if2),ytab(if2)
     $               ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                    ,F2Vave(if2),F2Eave(if2)



                  endif
C Search for same Q2:
                  do if22=if2+1,NF2TOT
                     if (.not.F2FLAG(if22)
     $                    .and.REACTAB(if22).eq.RLIST(irea)) then
                        if (Q2loc.eq.Q2TAB(if22)) then
                           nq2loc = nq2loc + 1
                           F2FLAG(if22) = .true.



                           if (mod(IOUT,10).eq.1) then
                              j = 0
                              do iexp=1,NMEAS
                                 if (RLIST(irea).eq.IREACT(iexp)) then
                                    j = j + 1
                                    f2tmp(j)  = f2tab(if22,iexp)
                                    f2tmpe(j) = f2etot(if22,iexp)
                                 endif
                              enddo
                              write (51,17) xtab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)
C All Q2 info:
                              write (50,117) q2loc,xtab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)
                           elseif (mod(IOUT,10).eq.2) then
C Float each point according to syst:
                              j = 0
                              do iexp=1,NMEAS
                                 if (RLIST(irea).eq.IREACT(iexp)) then
                                    j = j + 1
                                    f2tmp(j)  = f2tab(if22,iexp)
                                    f2tmpe(j) = f2etotn(if22,iexp)
                                    do isys=1,NSYSTOT
                                       f2tmp(j) = f2tmp(j) + 
     $                                      SYSTAB(isys,if22,iexp)
     $                                      *SYSSH(isys)
                                    enddo
                                 endif
                              enddo
                              write (51,17) xtab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)
C All Q2 info:
                              write (50,117) q2loc,xtab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)

C
C Add y ....
C
                           elseif (mod(IOUT,10).eq.3) then

                              j = 0
                              do iexp=1,NMEAS
                                 if (RLIST(irea).eq.IREACT(iexp)) then
                                    j = j + 1
                                    f2tmp(j)  = f2tab(if22,iexp)
                                    f2tmpe(j) = f2etot(if22,iexp)
                                 endif
                              enddo
                              write (51,17) xtab(if22),ytab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)
C All Q2 info:
                              write (50,117) q2loc,xtab(if22),ytab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)


C
C Add y ...
C
                           elseif (mod(IOUT,10).eq.4) then
C Float each point according to syst:
                              j = 0
                              do iexp=1,NMEAS
                                 if (RLIST(irea).eq.IREACT(iexp)) then
                                    j = j + 1
                                    f2tmp(j)  = f2tab(if22,iexp)
                                    f2tmpe(j) = f2etotn(if22,iexp)
                                    do isys=1,NSYSTOT
                                       f2tmp(j) = f2tmp(j) + 
     $                                      SYSTAB(isys,if22,iexp)
     $                                      *SYSSH(isys)
                                    enddo
                                 endif
                              enddo
                              write (51,17) xtab(if22),ytab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)
C All Q2 info:
                              write (50,117) q2loc,xtab(if22),ytab(if22)
     $                             ,(f2tmp(iexp),f2tmpe(iexp),iexp=1,j)
     $                             ,F2Vave(if22),F2Eave(if22)


                           endif
                        endif
                     endif
                  enddo
C Write Q2 info:
                  write (52,'(F8.2,i4)') q2loc,nq2loc
                  close (51)
               endif
 99         enddo
            close (50)
            close (52)
         enddo
      endif

C------------------------------------------------------------
      end
