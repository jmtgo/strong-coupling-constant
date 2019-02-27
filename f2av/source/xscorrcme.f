      Subroutine XSCORRCME(Reaction,q2,x,CMEIN,CMEOUT,XS)
C
C Created 7 Mar 05 by SG.
C Correct cross section for difference in Ceneter of mass energy:
C
C INPUT:
C   Reaction 
C   Q2,X
C   CMEIN
C   COMEPX
C INPUT/OUTPUT:
C   XS
C 
C Correction method:
C   XSout = XS + (XSpred(CMEOUT)-XSperd(CMEIN))
C
      implicit none
      include '../include/f2.inc'
      integer Reaction
      real q2,x,CMEIN,CMEOUT,XS
C
      Real F2a,F2b

C
C -- H1 parametrizations of NC/CC cross sections:
C
      real signc,sigcc,signc2,FitFract

C--------------
      if (CMEOUT.ne.CMEIN) then
         if (Reaction.eq.515) then
C NC e- p

C 31 may 05 SG: no correction for low Q2 at the moment:
            if (Q2.gt.4.) then
               F2a = signc2(x,q2,-1.,CMEIN)
               F2b = signc2(x,q2,-1.,CMEOUT)
            else
               F2a = 1.
               F2b = 1.
            endif
         elseif (Reaction.eq.615) then
C NC e+ p
            if (LFractal) then
               F2a = FitFract(x,q2,q2/(CMEIN*x))
               F2b = FitFract(x,q2,q2/(CMEOUT*x))
            else
               F2a = signc2(x,q2,+1.,CMEIN)
               F2b = signc2(x,q2,+1.,CMEOUT)
            endif
         elseif (Reaction.eq.3515) then
C CC e- p
            F2a =  sigcc(x,q2,-1.,CMEIN)
            F2b =  sigcc(x,q2,-1.,CMEOUT)
            
         elseif (Reaction.eq.3615) then
C CC e+ p
            F2a =  sigcc(x,q2,+1.,CMEIN)
            F2b =  sigcc(x,q2,+1.,CMEOUT)
            
         else
            print *,'UNKNOWN REACTION=',REACTION
            print *,'NO CoME correction'
            F2a = 1.
            F2b = 1.
         endif
         if (abs(f2a-f2b)/XS.gt.0.002) then
            print 77,q2,
     $           q2/x/cmein,q2/x/CMEOUT,f2a,f2b,100.*(f2b-f2a)/XS
 77         format('CME >0.2%',F10.2,2F6.2,2F7.3,F7.1)
         endif

         if (Reaction.eq.515 .or. Reaction.eq.615) then
            
            if (LOGAVE) then
               XS = log ( exp(XS) - f2a + f2b)
            else
               XS = XS - f2a + f2b
            endif
         elseif (Reaction.eq.3515 .or. Reaction.eq.3615) then
            if (LOGAVE) then
               XS = XS + log(F2b/F2a)
            else
               XS = XS * (F2b/F2a)
            endif
         endif
      endif
      end
