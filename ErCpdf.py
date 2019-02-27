'''This module evaluates the Cpdf matrix as a percentage error (rather than as an absolute error,
it takes in values of shelves generated in Hessian.py or NNPDF.py This matrix for each pdf set can be used 
to reverse engineer the problem to find the error matrix for the non-principal as values'''

import numpy as np
import shelve


Er = shelve.open('ErPdf')
yNN = shelve.open('yNN')


def run_scale_variation():
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
            if a == 1 and b == 1:
                # Do nothing in this case
                print a,b
            else:
                print a,b
                theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_('+str(a)+',' +str(b)+ ')')
                theory_z3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet_('+str(a)+',' +str(b)+ ')')

                Cpd2 = shelve.open('Cpdf2_('+str(a)+',' +str(b)+ ')') #will need to specify the PDF later in the code
                Cpd3 = shelve.open('Cpdf3_('+str(a)+',' +str(b)+ ')') #(1,1) (1,0.5) (1,2) (0.5,0.5) (0.5,1) (0.5,2) (2,1) (2,2) (2,0.5)
                

                #Er['CT10nlo2_('+str(a)+',' +str(b)+ ')'] = Cov_100(theory_z2, new_pdf_name, numberjet_0th_value_best_a_s, y2, Cpd2)
                #Er['CT10nlo3_('+str(a)+',' +str(b)+ ')'] = Cov_100(theory_z3, new_pdf_name, numberjet_0th_value_best_a_s, y3, Cpd3)
                for i in range(118, 122):
                    values = Cov_allas(theory_z3, 'CT14_nlo_0'+str(i)+'', y3, 'CT14nlo3_('+str(a)+',' +str(b)+ ')')
                    Cpd3['CT14nlo0' + str(i) + ''] = values


                for i in range(118,119):
                    values = Cov_allas(theory_z2, 'CT14_nlo_0'+str(i)+'', y2, 'CT14nlo2_('+str(a)+',' +str(b)+ ')')
                    Cpd2['CT14nlo0' + str(i) + ''] = values

                theory_z2.close()
                theory_z3.close()
                #theory_z4.close()
                Cpd2.close()
                Cpd3.close()
                #Cpd4.close()

def si(jet, label,y): #gives the fastnlo cross section bin data.
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet[label])[1:y+1] #for particular pdf parameter change k, is a vector length y


def CovNN_100(Cplabel, y, Cpdjet, yNN):
    #finds what the error is by "reverse percentaging" of the actual measurment
    X = np.zeros((y,y), dtype= np.float64)
    for i in np.arange(y):
        for j in np.arange(y):
            X[i,j] = (Cpdjet[Cplabel][i,j]) *100 / (yNN[i]*yNN[j]) 
    return X

def Cov_100(jet, Cplabel, label, y, Cpdjet):
    #finds what the error is by "reverse percentaging" of the actual measurment
    X = np.zeros((y,y), dtype= np.float64)
    for i in np.arange(y):
        for j in np.arange(y):
            X[i,j] = (Cpdjet[Cplabel][i,j]) *100 / (si(jet,label,y)[i]*si(jet,label,y)[j])
    return X

def  Cov_allas(jet, label, y, Erpdf): #label is the label of the as pdf, Erpdf is the label of the Erpdf shelve.
    X = np.zeros((y,y), dtype= np.float64)
    for i in np.arange(y):
        for j in np.arange(y):
            X[i,j] = Er[Erpdf][i,j] /100 * (si(jet,label,y)[i]*si(jet,label,y)[j])
    return X 

if __name__ =="__main__":

    y2 = 13
    y3 = 13
    y4 = 12  
    x = np.zeros(y4)
    z = np.zeros(y4)



    '''REMEMBER TO CHANGE THE EPDF SHELF ABOVE'''

    run_scale_variation()

    #Er['MSTW082_('+str(a)+',' +str(b)+ ')'] = Cov_100(t2, 'MSTW08cl', 'MSTW0868cl0', y2, Cpd2)
    #Er['MSTW083_('+str(a)+',' +str(b)+ ')'] = Cov_100(t3, 'MSTW08cl', 'MSTW0868cl0', y3, Cpd3)
    #for i in range(y4):

    #    x[i] = (Cpd2['NN23_0120'])[i,i]
    #print np.sqrt(x)

    #for i in range(y4):

   #     z[i] = (Cpd2['MSTW08clo012'])[i,i]
   # print np.sqrt(z)
    
    #for i in range(0, 22):
    #    values = Cov_allas(t3, 'MSTW08_000'+str(i)+'', y3, 'MSTW083_('+str(a)+',' +str(b)+ ')')
    #    Cpd3['MSTW08clo0' + str(i) + ''] = values


 



Er.close()
yNN.close()
