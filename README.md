# strong-coupling-constant

The purpose of this project is to establish a value for the strong coupling constant and the associated error.
Cross section values are caluculated for different pdf sets from LHAPDF using fastnlo.

To calculate the error we will need to calculate the covariance matrix with Cpdf, Cstat, Ctheory components.
Cstat is calculated from the errors generated in the HEP data set. 

Cpdf will be calculated using data from pdf sets listed below:

These pdfs will use Hessian.py:
- CT10nlo, CT10nlo_as_0i
- CT14nlo, CT14nlo_as_0i
- MSTW2008nlo68cl, MSTW2008nlo68cl_asmzrange
- MMHT2014nlo68cl, MMHT2014nlo68cl_asmzrange

These pdfs will use a statistical approach:
- NNPDF23_nlo_as_0i
- NNPDF30_nlo_as_0i

twojet, threjet and fourjet are shelves with different keys with the following structure:
[0112] = CT10 pdf set with as = 0112
[CT10nloij] = CT10 pdf set 2nd approximation for the pdf up to 68% confidence level 'ij' number.
[CT14nlo_as_0ijk] =CT14 pdf set with as = ijk
[CT14nloij] = ct14 PDF SET, 2nd approximation for the pdf up to 68% confidence level 'ij' number.


