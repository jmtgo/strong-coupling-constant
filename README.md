# strong-coupling-constant

The purpose of this project is to establish a value for the strong coupling constant and the associated error.
Cross section values are caluculated for different pdf sets from LHAPDF using fastnlo.

To calculate the error we will need to calculate the covariance matrix with Cpdf, Cstat, Ctheory components.
Cstat is calculated from the errors generated in the HEP data set. 

Cpdf will be calculated using data from pdf sets listed in linux-shelves read me.

twojet, threjet and fourjet are shelves with different keys, they are stored in the repositary linux-shelves.

for the neural network pdf's:
- [NNPDFij_0k_l] = ij = nn type, k = value of as, l = index of the pdf set.

The Cpdf's were calculated using Hessian and NNPDF (for the NNPDF pdf sets), and were stored in Cpdfi, where i is the jet number: The keys of which are:
- [MMHT14cl] = MMHT14 for the principle value of as.
- [MSTW08cl] = MSTW08 for the principle value of as. 
- [CT10nlo] = CT10nlo for the principle value of as.
- [CT14nlo] = CT14nlo for the principle value of as.
- [NN23_0ii] = NNPDF23 for the value of as = 0ii. (114-124)
- [NN30_0ii] = NNPDF30 for the value of as = 0ii. (115,117,118,119,121).

The Cpdf's were converted into percentage error matrices and stored in ErPdf:
- [CT10nloi] = jet i principle value of as.
- [CT14nloi] = jet i principle value of as.
- [MMHT14i] = jet i principle value of as.
- [MSTW08i] = jet i principle value of as.

These were then used to 'reverse' calculate for the non principle as values and were stored in Cpdfi also, with keys:
- [MMHT14clo0ii] = MMHT for the ii value of as. (0.108-0.128) (1,20) 0TH IS PRINCIPLE
- [MSTW08cl0ii] = MSTW for the ii value of as. (0.107-0.128) (1,21), 0TH IS PRINCIPLE
- [CT10nlo0ii] = CT10 for the ii value of as. (0.112-0.128) (112,127)
- [CT14nlo0ii] = CT14 for the ii value of as. (0.111-0.124) (111,124)

The average values for some as values were calculated for NNPDF, and stored in yNN:
-[NN30_0iii_jk] = for NNPDF30 as = iii (118 so far), k = jet number.




