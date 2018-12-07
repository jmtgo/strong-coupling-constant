# strong-coupling-constant

The purpose of this project is to establish a value for the strong coupling constant and the associated error.
Cross section values are caluculated for different pdf sets from LHAPDF using fastnlo.

To calculate the error we will need to calculate the covariance matrix with Cpdf, Cstat, Ctheory components.
Cstat is calculated from the errors generated in the HEP data set. 

Cpdf will be calculated using data from pdf sets listed in linux-shelves read me.

For different factorisation and renormalisation values there is an individual shelve, with format twojet_(0.5,0.5), for all combinations of 0.5, 1, 2, except (1,1) which is just twojet as this is the best approximation. The keys within each shelve is the same. There are eight scales in total (Except the (1,1) scale). 

twojet, threjet and fourjet are shelves with different keys with the following structure, they are stored in the repositary linux-shelves:
- [0112] = CT10 pdf set with as = 0112
- [CT10nloij] = CT10 pdf set 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [CT14_nlo_0ijk] =CT14 pdf set with as = ijk
- [CT14nloij] = ct14 PDF SET, 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [MMHT14_000ij] = MMHT PDF SET with as = ij
- [MMHT1468cli] = MMHT PDF SET, 2nd approximation for the pdf up to 68% confidence level 'ij' number.
- [MSTW08_000ij] = MSTW PDF SET with as = ij
- [MSTW0868cli] = MSTW PDF SET, 2nd approximation for th pdf up to 68% confidence level 'ij' number.


for the neural network pdf's:
- [NNPDFij_0k_l] = ij = nn type, k = value of as, l = index of the pdf set.

for the neural network pdf's:
- [NNPDFij_0k_l] = ij = nn type, k = value of as, l = index of the pdf set.

The Cpdf's were calculated using Hessian and NNPDF (for the NNPDF pdf sets), and were stored in Cpdfi, where i is the jet number. The keys of which are for (1,1):
- [MMHT14cl] = MMHT14 for the principle value of as.
- [MSTW08cl] = MSTW08 for the principle value of as. 
- [CT10nlo] = CT10nlo for the principle value of as.
- [CT14nlo] = CT14nlo for the principle value of as.
- [NN23_0ii] = NNPDF23 for the value of as = 0ii. (114-124)
- [NN30_0ii] = NNPDF30 for the value of as = 0ii. (115,117,118,119,121).

The Cpdf's were calculated using Hessian and NNPDF (for the NNPDF pdf sets), and were stored in Cpdfi_(a,b), where i is the jet number, (a,b) is the scale, with the same keys as for Cpdfi with scale (1,1).

The Cpdf's were converted into percentage error matrices and stored in ErPdf, for scale (1,1):
- [CT10nloi] = jet i principle value of as.
- [CT14nloi] = jet i principle value of as.
- [MMHT14i] = jet i principle value of as.
- [MSTW08i] = jet i principle value of as.

For scales of all combinations besides (1,1), the ErPdf: 

- [MSTW08i_(0.5,0.5)] = jet i principle value of as scale (0.5,0.5)

These were then used to 'reverse' calculate for the non principle as values and were stored in Cpdfi also, with keys:
- [MMHT14cl0ii] = MMHT for the ii value of as. (0.108-0.128) (1,20) 0TH IS PRINCIPLE
- [MSTW08clo0ii] = MSTW for the ii value of as. (0.107-0.128) (1,21), 0TH IS PRINCIPLE
- [CT10nlo0ii] = CT10 for the ii value of as. (0.112-0.127) (112,127)
- [CT14nlo0ii] = CT14 for the ii value of as. (0.111-0.123) (111,123)

The average values for some as values were calculated for NNPDF, and stored in yNN (for (1,1)):
- [NN30_0iii_jk] = for NNPDF30 as = iii (all as not jet 4), k = jet number.
- [NNPDF23_0iii_jk] = for NNPDF23 as=iii, k = jet number. (not done jet 4, done all as so far)




