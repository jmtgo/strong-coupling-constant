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
- Note the confidence level pdf sets come in pair, where the 0th set is the best estimate and i, i+1 are in pairs.

These pdfs will use NNPDF.py:
- NNPDF23_nlo_as_0i
- NNPDF30_nlo_as_0i

twojet, threjet and fourjet are shelves with different keys with the following structure, they are stored in the repisitary linux-shelves:
- [0112] = CT10 pdf set with as = 0112
- [CT10nloij] = CT10 pdf set 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [CT14nlo_as_0ijk] =CT14 pdf set with as = ijk
- [CT14nloij] = ct14 PDF SET, 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [MMHT14_000ij] = MMHT PDF SET with as = ij
- [MMHT1468cli] = MMHT PDF SET, 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [MSTW08_000ij] = MSTW PDF SET with as = ij
- [MSTW0868cli] = MSTW PDF SET, 2nd approximation for th pdf up to 68% confidence level 'ij' number.


for the neural network pdf's:
- [NNPDFij_0k_l] = ij = nn type, k = value of as, l = index of the pdf set.

The Cpdf's were calculated using Hessian and NNPDF (for the NNPDF pdf sets), the keys of which are:
- [MMHT14cl] = MMHT14 for the principle value of as.
- [MSTW08cl] = MSTW08 for the principle value of as.
- [CT10nlo] = CT10nlo for the principle value of as.
- [CT14nlo] = CT14nlo for the principle value of as.
- [NN23_0ii] = NNPDF23 for the value of as = 0ii. (114-124)
- [NN30_0ii] = NNPDF30 for the value of as = 0ii. (115,117,118,119,121).



