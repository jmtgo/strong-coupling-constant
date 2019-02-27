      Program F2AVE
C-----------------------------------------------------------
C
C Created 4 June 04 by SG.
C
C usage: f2ave.exe < steer
C
C Read in F2(x,Q), Sigma(F2), dF2/dSyst and Sigma(Syst) from
C different experiments and calculate average/covariance matrix 
C
C-----------------------------------------------------------
      implicit none
      include '../include/f2.inc'
      
C--------------------------------------------------------------------
      print *,'--------------------------------'
     $     //'---------------------------------------------------'
      print *,'    Initiating Linear Average With Systematic '//
     $     'UncertaInTies  program'
      print *,' '
      print *,' Version of Thu Sep 10 2009 '
      print *,' '
      print *,'---------------------------------------------------'
     $     //'--------------------------------'

      Call readf2         ! Read x-section tables
      if (NSIM.gt.0) then
         Call simf2
      endif

      Call fitf2          ! Perform F2 averaging
      Call outf2          ! Write additional averaging info

C--------------------------------------------------------------------
      end
