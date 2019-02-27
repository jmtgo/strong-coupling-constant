*
* ..File: cq1f2lh.f
*
* ..The tabulated NLO light-quark coefficient functions for the heavy 
*    quark contributions to  F_L  and  F_T = F_2 - F_L.  This is an 
*    improved version of the program presented in
*
*    Riemersma, Smith and van Neerven, PL B347 (1995) 143  (RSN).
*
* ..The functions  clq1, ctq1  [dlq1, dtq1]  give the scale-independent
*    parts proportional to  e_h^2 [e_l^2], respectively;  clq1l  and 
*    ctq1l  the corresponding coefficients of  ln (mu_f^2/m_h^2).  
*    The additional mass-factorization for Q^2 -> 0 in the e_l^2 part
*    is not included.  The prefactors are  a_s^2 = (alpha_s/4 pi)^2,
*    unlike in RSN.
*    
*    eta = (W^2 - 4.d0*m2)/4.d0/m^2  where W is the CM energy of the
*    gamma*parton system. xi = Q^2/m^2
*
* =====================================================================
*
*
      function clq1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
*
      clq1 = cf*tf * ( beta**3 * asymp_l(xi) + hlq_h(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      clq1 = 16.d0 * pi * xi * clq1
*
      return
      end
*
* ---------------------------------------------------------------------
*
      function ctq1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
*
      ctq1 = cf*tf * ( beta**3 * asymp_t(xi) + htq_h(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      ctq1 = 16.d0 * pi * xi * ctq1
*
      return
      end
*
* ---------------------------------------------------------------------
*
      function clq1l (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
*
      clq1l = cf*tf * ( beta**3 * asympb_l(xi) + hblq_h(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      clq1l = 16.d0 * pi * xi * clq1l
*
      return
      end
*
* ---------------------------------------------------------------------
*
      function ctq1l (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
*
      ctq1l = cf*tf * ( beta**3 * asympb_t(xi) + hbtq_h(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      ctq1l = 16.d0 * pi * xi * ctq1l
*
      return
      end
*
* ---------------------------------------------------------------------
*
      function dlq1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      dlq1 = cf*tf * hlq_l(eta,xi) 

C ... Comment out the following line (from AV) to get the original NPB function
C      dlq1 = 16.d0 * pi* xi * cf*tf * hlq_l(eta,xi) 
*
      return
      end
*
* ---------------------------------------------------------------------
*
      function dtq1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      dtq1 = cf*tf * htq_l(eta,xi) 

C ... Comment out the following line (from AV) to get the original NPB function
C      dtq1 = 16.d0 * pi* xi * cf*tf * htq_l(eta,xi) 
*
      return
      end
*
*
* =================================================================av==
