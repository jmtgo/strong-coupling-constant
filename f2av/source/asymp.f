*
*    File: asymp.f (from SR)
*
*    These are the functions that give the asymptotic dependence of the
*    coefficient functions with the appropriate factors.
*
*    xi = Q^2/m2
*
* ======================================================================
*
*    longitudinal: equation (19) in PLB347 (1995) 143 - 151
*
      real*8 function asymp_l(xi)
      implicit none
      real*8 xi, pi, term1, fii, fjj
      parameter (pi = 3.14159265359d0)
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      asymp_l = 1.d0/6.d0/pi*(4.d0/xi - 4.d0/3.d0*term1 
     #     + (1.d0 - 2.d0/xi - 1.d0/6.d0*term1)*fjj(xi)
     #     - (3.d0/xi + 0.25d0*term1)*fii(xi))
      return
      end
*
c transverse: equation (20) 
*
      real*8 function asymp_t(xi)
      implicit none
      real*8 xi, pi, term1, fii, fjj
      parameter (pi = 3.14159265359d0)
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      asymp_t = 1.d0/6.d0/pi*(-2.d0/3.d0/xi + 4.d0/3.d0*term1
     #     + (7.d0/6.d0 + 1.d0/3.d0/xi + 1.d0/6.d0*term1)*fjj(xi)
     #     + (1.d0 + 2.d0/xi + 0.25d0*term1)*fii(xi))
      return
      end
*
c longitudinal mass factorization: equation (21) 
*
      real*8 function asympb_l(xi)
      implicit none
      real*8 xi, pi, term1, fjj
      parameter (pi = 3.14159265359d0)
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      asympb_l = 1.d0/6.d0/pi*(-6.d0/xi + 0.5d0*term1
     #     + (3.d0/xi + 0.25d0*term1)*fjj(xi))
      return
      end
*
c transverse mass factorization: equation (22) 
*
      real*8 function asympb_t(xi)
      implicit none
      real*8 xi, pi, term1, fjj
      parameter (pi = 3.14159265359d0)
      term1 = 1.d0/(1.d0 + 0.25d0*xi)
      asympb_t = 1.d0/6.d0/pi*(4.d0/xi - 0.5d0*term1
     #     - (1.d0 + 2.d0/xi + 0.25d0*term1)*fjj(xi))
      return
      end
*
c equation (23) 
*
      real*8 function fjj(xi)
      implicit none
      real*8 pi, xi, term1, term2
      parameter (pi = 3.14159265359d0)
      term1 = dsqrt(xi)
      term2 = dsqrt(4.d0 + xi)
      fjj = 4.d0/term1/term2*dlog((term2 + term1)/(term2 - term1))
      return
      end
*
c equation (24) 
*
      real*8 function fii(xi)
      implicit none
      real*8 pi, term1, term2, xi, dilog
      parameter (pi = 3.14159265359d0)
      term1 = dsqrt(xi)
      term2 = dsqrt(4.d0 + xi)
      fii = 4.d0/term1/term2*(-pi*pi/6.d0 
     #      - 0.5d0*(dlog((term2 + term1)/(term2 - term1)))**2
     #      + (dlog(0.5d0*(1.d0 - term1/term2)))**2 
     #      + 2.d0*dilog(0.5d0*(1.d0 - term1/term2)))
      return
      end
*
*=======================================================================
