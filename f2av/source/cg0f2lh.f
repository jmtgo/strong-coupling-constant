*
*    File cg0f2lh.f (from SR, normalization changed and renamed by AV)
*
*    This gives the transverse Born coefficient function as shown in 
*    fig 6a of NPB392 (1993) 162 - 229.  For QCD take tf = 1d0/2d0, 
*    for QED take tf = 1d0. 
*
*    eta = (W^2 - 4d0*m2)/4d0/m^2  where W is the CM energy of the
*    gamma*parton system. xi = Q^2/m^2
*
* ======================================================================
*
      function ctg0(eta,xi)
*
      implicit none
      real*8 eta, xi, pi, ca, cf, tf, ctg0
      common/group/ca, cf, tf
      parameter(pi = 3.14159265359d0)
*
      ctg0 = 0.5d0*pi*tf*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #       (-2.d0*((1.d0 + eta - 0.25d0*xi)**2 + eta + 1.d0)*
     #       sqrt(eta/(1.d0 + eta)) + (2.d0*(1.d0 + eta)**2 +
     #       0.125d0*xi**2 + 2.d0*eta + 1.d0)*
     #       log((sqrt(1.d0 + eta) + sqrt(eta))/
     #           (sqrt(1.d0 + eta) - sqrt(eta))))

C ... Comment out the following line (from AV) to get the original NPB function
C      ctg0 = xi/pi * ctg0 
*
      return
      end
*
c Longitudinal coefficient function, fig 6b of NPB392(1993) 162 - 229.
*
      function clg0(eta,xi)
*
      implicit none
      real*8 eta, xi, pi, tf, ca, cf, clg0
      common/group/ca, cf, tf
      parameter(pi = 3.14159265359d0)
*
      clg0 = 0.5d0*pi*tf*xi*(1.d0 + eta + 0.25d0*xi)**(-3)*
     #       (2.d0*sqrt(eta*(1.d0 + eta)) -
     #       log((sqrt(1.d0 + eta) + sqrt(eta))/
     #           (sqrt(1.d0 + eta) - sqrt(eta))))

C ... Comment out the following line (from AV) to get the original NPB function
C      clg0 = xi/pi * clg0 
*
      return
      end
*
* ======================================================================
