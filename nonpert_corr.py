'''This module evaluates the covariance matrix for the non-perturbative 
corrections between hadron and parton level. Also includes the corrections to the predictions'''

import shelve
import numpy as np
import pandas as pd


jet2 = pd.read_csv('HEPData-ins1514251-v2-Table_66.csv')
jet3 = pd.read_csv('HEPData-ins1514251-v2-Table_67.csv')


#nocorr2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_('+str(a)+',' +str(b)+ ')')
#nocorr3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet_('+str(a)+',' +str(b)+ ')')
#for the non neural network pdfs we use the cross sections straight from fastnlo, for nn we use yNN shelves.



#Cpd2 = shelve.open('Cpdf2__('+str(a)+',' +str(b)+ ')') #only done for principle factorisation values. 
#Cpd3 = shelve.open('Cpdf3__('+str(a)+',' +str(b)+ ')') #(1,1) (1,0.5) (1,2) (0.5,0.5) (0.5,1) (0.5,2) (2,1) (2,2) (2,0.5)


def run_scale_variation(pdf, start,stop, step):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
            #the yNN values are presliced.
            nocorr = shelve.open('yNN_('+str(a)+',' +str(b)+ ')')

            theory_cross2_full = shelve.open('pc2__('+str(a)+',' +str(b)+ ')')
            theory_cross3_full = shelve.open('pc3__('+str(a)+',' +str(b)+ ')')

            pred_correc(pdf, start, stop, step, '_j2', theory_cross2_full, corr2, nocorr)
            pred_correc(pdf, start, stop, step, '_j3', theory_cross3_full, corr3, nocorr)

            nocorr.close()
            theory_cross2_full.close()
            theory_cross3_full.close()



def si(jet_number_pdf, pdf,y,i):
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet_number_pdf[pdf+str(i) + ''])[1:y+1] #for particular pdf parameter change k, is a vector length y


def Cov_non_pert_correc(y2, y3, ecorr2, ecorr3, pdf, start, stop, step):
    a=1
    b=1
    for i in range(start, stop, step):
        ypred = shelve.open('ypred_('+str(a)+','+str(b)+')')
        X = np.zeros(((y2+y3)*2,(y2+y3)*2), dtype=np.float64) #not fully correlated, try both diagonals and this way. 
        X = np.diag(full_pred_error( ecorr2, ecorr3)**2) * ypred[pdf+str(i)]
        print np.shape(X)

    return X

    ypred.close()

def pred_correc(pdf, start, stop, step, j, pc, corr, nocorr):
    
    for i in range(start, stop, step):

        values = corr * nocorr[pdf+str(i) + j]
        print values
        pc[pdf +str(i)+ ''] = values
    print len(list(pc.keys()))
    

def full_pred_error( ecorr2, ecorr3):
    A = np.append([ecorr2],[ecorr3])
    B =  np.append(A, [ecorr2] )
    C = np.append(B, [ecorr3] )
    return C
    




if __name__ =="__main__":
    y2 = 13
    y3 = 13
    y4 = 12   
    ecorr2 = np.array(jet2['syst +'])
    ecorr3 = np.array(jet3['syst +'])
    corr2 = np.array(jet2['NP_corrections_Zll [none]'])
    corr3 = np.array(jet3['NP_corrections_Zll [none]'])


    pdf = 'NNPDF30_0'
    start = 118
    stop = 119
    step = 1
    #run_scale_variation(pdf, start,stop, step)
    Cov_non_pert_correc(y2, y3, ecorr2, ecorr3, pdf, start, stop, step)
    
    #Cpd2['pert_corr'] = Cov_non_pert_correc(y2, corr2)
    #Cpd3['pert_corr'] = Cov_non_pert_correc(y3,corr3)
   