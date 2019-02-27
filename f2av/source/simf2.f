      Subroutine SimF2
*
*
* Created on 29 Mar 2007 by SG
*
* Simulate X-section based on read tables
*
*
      implicit none
      include '../include/f2.inc'
      integer if2,iexp,isys

      real systsh(NSYSTMAX)
      real f2sh(NF2MAX),f2orig,f2new,scale
      integer ntotin,nto2in,ioffs
C------------------------------------------------
      print *,'SIMULATE INPUT TABLES'



      if (NMEAS.gt.1) then
         print *,'Currently implemented for 1 input dataset only'
         print *,'Aborting'
         stop
      endif

C initialize rnd numbers:
      call rmarin(ISEED,ntotin,nto2in)


      do iexp=NMEAS+1,NMEAS+NSIM

         ioffs = (iexp-NMEAS)*10

C random shift of systematics:
         call RNORML(systsh(isys),NSYSTOT)
C random shift of XSection:
         call RNORML(f2sh,NF2TOT)
C Apply
         do if2=1,NF2TOT
            f2orig = F2TAB(if2,1)
            f2new = f2orig + F2ETAB(if2,1)*f2sh(if2)
            do isys=1,NSYSTOT
               f2new = f2new + SYSTAB(isys,if2,1)*systsh(isys)
            enddo
            scale = f2new/f2orig
C Store:
            F2TAB(if2,iexp)  = f2new
            F2ETAB(if2,iexp) = F2ETAB(if2,1)*Scale
            F2ETABOrig(if2,iexp) =   F2ETAB(if2,1)*Scale
            do isys=1,NSYSTOT
               SYSTAB(isys+ioffs,if2,iexp) = SYSTAB(isys,if2,1)*scale
               SYSTABorig(isys+ioffs,if2,iexp)
     $              = SYSTAB(isys,if2,1)*scale
            enddo
         enddo         
      enddo
      NSYSTOT = NSYSTOT + NSIM*10
      NMEAS = NMEAS + NSIM
      end
