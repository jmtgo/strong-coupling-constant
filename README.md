# strong-coupling-constant

The purpose of this project is to establish a value for the strong coupling constant and the associated error.
Cross section values are caluculated for different pdf sets from LHAPDF using fastnlo.

To calculate the error we will need to calculate the covariance matrix with Cpdf, Cstat, Ctheory components.
Cstat is calculated from the errors generated in the HEP data set. 

Cpdf will be calculated using data from pdf sets listed in linux-shelves read me.

twojet, threjet and fourjet are shelves with different keys, they are stored in the repositary linux-shelves.

for the neural network pdf's:
- [NNPDFij_0k_l] = ij = nn type, k = value of as, l = index of the pdf set.

The Cpdf's were calculated using Hessian and NNPDF (for the NNPDF pdf sets), the keys of which are:
- [MMHT14cl] = MMHT14 for the principle value of as.
- [MSTW08cl] = MSTW08 for the principle value of as.
- [CT10nlo] = CT10nlo for the principle value of as.
- [CT14nlo] = CT14nlo for the principle value of as.
- [NN23_0ii] = NNPDF23 for the value of as = 0ii. (114-124)
- [NN30_0ii] = NNPDF30 for the value of as = 0ii. (115,117,118,119,121).

The Cpdf's were converted into percentage error matrices and stored in ErPdf:
- [CT10nloi] = jet i principle value of as.
- [CT14nloi[ = jet i principle value of as.




