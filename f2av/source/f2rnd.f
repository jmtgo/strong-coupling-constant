      Subroutine F2RND(F2Ran,Sysran,imode)
C
C Generate random F2 and systematics
C
      implicit none
      include '../include/f2.inc'
      real*8 F2Ran(NF2MAX),Sysran(NSYSTMAX)
      integer imode
      integer i,j
      real RNDM
C----------------------------------------------
      

      if (imode.ne.0) then
         do i=1,NF2TOT
            F2RAN(i) = 1.5 + (RNDM()-0.5)*0.2
         enddo
         do i=1,NSYSTOT
            SYSRAN(i) = (RNDM()-0.5)*2
         enddo
         
      else
C Get first suitable F2:
         do i=1,NF2TOT
            do j=1,NMEAS
               if (F2TAB(i,j).ne.0) then
                  F2RAN(i) = F2TAB(i,j)
                  goto 17
               endif
            enddo
 17         continue
         enddo
         
         do i=1,NSYSTOT
            SYSRAN(i) = 0.0 ! (RNDM()-0.5)*0.01
         enddo

      endif



      end
