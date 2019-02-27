*
* File: thresh.f (from SR)
*
* These are the functions that give the threshold dependence of the 
* coefficient functions with the appropriate factors.
* eta = (W^2 - 4d0*m2)/4d0/m^2  where W is the CM energy of the gamma* 
* parton system. xi = Q^2/m^2
*
c Longitudinal CF group structure: equation (13) in PLB347 (1995) 143
*
      double precision function thrshf_l(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrshf_l = 1.d0/6.d0/pi*xi*term1**3*beta*beta*pi*pi/2.d0
      return
      end
*
c Transverse CF group structure: equation (14) 
*
      double precision function thrshf_t(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrshf_t = 0.25d0/pi*term1*pi*pi/2.d0
      return
      end
*
c Longitudinal CA group structure: equation (15)
*
      double precision function thrsha_l(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrsha_l = 1.d0/6.d0/pi*xi*term1**3*beta**2*
     #     (beta*(dlog(8.d0*beta*beta))**2
     #     - 5.d0*beta*dlog(8.d0*beta*beta) - 0.25d0*pi*pi)
      return
      end
*
c Transverse CA group structure: equation (16) 
*
      double precision function thrsha_t(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrsha_t = 0.25d0/pi*term1*(beta*(dlog(8.d0*beta*beta))**2
     #     - 5.d0*beta*dlog(8.d0*beta*beta) - 0.25d0*pi*pi)
      return
      end
*
c Longitudinal CA group structure for the mass factorization piece: 
c equation (17) in PLB347 (1995) 143 - 151
*
      double precision function thrshb_l(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrshb_l = 1.d0/6.d0/pi*xi*term1**3*beta**3*
     #     (-dlog(4.d0*beta*beta))
      return
      end
*
c Transverse CA group structure for the mass factorization piece: 
c equation (18) in PLB347 (1995) 143 - 151
*
      double precision function thrshb_t(eta,xi)
      implicit none
      double precision pi, eta, xi, beta, term1
      parameter (pi = 3.14159265359d0)
      beta = dsqrt(eta/(1.d0 + eta))
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      thrshb_t = 0.25d0/pi*term1*beta*(-dlog(4.d0*beta*beta))
      return
      end
*
*=======================================================================
