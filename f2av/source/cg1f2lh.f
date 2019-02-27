*
* ..File: cg1f2lh.f
*
* ..The tabulated NLO gluon coefficient functions for the heavy quark
*    contributions to  F_L  and  F_T = F_2 - F_L.  This is an improved 
*    version of the program presented in
*
*    Riemersma, Smith and van Neerven, PL B347 (1995) 143  (RSN).
*
* ..The functions  clg1 and ctg1  return the scale-independent parts, 
*     clg1l and ctg1l  the coefficients of  ln (mu_f^2/m_h^2).  The 
*     prefactors are simply  a_s^2 = (alpha_s/4 pi)^2,  unlike in RSN.
*    
*    eta = (W^2 - 4d0*m2)/4d0/m^2  where W is the CM energy of the
*    gamma*parton system. xi = Q^2/m^2
*
* =====================================================================
*
*     see equation (9) and (10) in PLB347 (1995) 143
*     
      function clg1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
      rho  = 1.d0/(1.d0 + eta)
*
      clg1 =  ca*tf * ( beta * asymp_l(xi) + rho * thrsha_l(eta,xi)
     +                + hlg_ca(eta,xi) )
     +      + cf*tf * ( rho * thrshf_l(eta,xi) + hlg_cf(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      clg1 = 16.d0 * pi * xi * clg1
*
      return
      end
*
* ---------------------------------------------------------------------
*
*     see equation (9) and (10) in PLB347 (1995) 143
*     
      function ctg1 (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
      rho  = 1.d0/(1.d0 + eta)
*
      ctg1 =  ca*tf * ( beta * asymp_t(xi) + rho * thrsha_t(eta,xi)
     +                + htg_ca(eta,xi) )
     +      + cf*tf * ( rho * thrshf_t(eta,xi) + htg_cf(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      ctg1 = 16.d0 * pi * xi * ctg1
*
      return
      end
*
* ---------------------------------------------------------------------
*
*     see equation (12) in PLB347 (1995) 143
*     
      function clg1l (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
      rho  = 1.d0/(1.d0 + eta)
*
      clg1l = ca*tf * ( beta* asympb_l(xi) + rho* thrshb_l(eta,xi)
     +                + hblg_ca(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      clg1l = 16.d0 * pi * xi * clg1l
*
      return
      end
*
* ---------------------------------------------------------------------
*
*     see equation (12) in PLB347 (1995) 143
*     
      function ctg1l (eta,xi)
*
      implicit real*8 (a -z)
      parameter (pi = 3.14159265359d0)
      common / group / ca, cf, tf
*
      beta = sqrt(eta/(1.d0 + eta))
      rho  = 1.d0/(1.d0 + eta)
*
      ctg1l = ca*tf * ( beta* asympb_t(xi) + rho* thrshb_t(eta,xi)
     +                + hbtg_ca(eta,xi) )

C ... Comment out the following line (from AV) to get the original NPB function
C      ctg1l = 16.d0 * pi * xi * ctg1l
*
      return
      end
*
* =================================================================av==
