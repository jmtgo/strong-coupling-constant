'''This module evaluates the Cpdf matrix for NN pdfs'''
import numpy as np
import shelve

theory_z2 = shelve.open('twojet')
theory_z3 = shelve.open('threejet')
theory_z4 = shelve.open('fourjet')
Cpd2 = shelve.open('Cpdf2')
Cpd3 = shelve.open('Cpdf3')
Cpd4 = shelve.open('Cpdf4')



def si(jet, label, y, k, i): #label = 'NNPDF23_0'
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    # k is the index of the alpha value, i is the index of the index of the pdf set
    return np.array(jet[label + str(k) + '_' + str(i)+''])[1:y+1] #for particular pdf parameter change k, is a vector length y

def Av(jet,label,y,k):
    Av = np.zeros((y), dtype=np.float64)
    for i in range(101):
        Av = si(jet, label, y, k, i) + Av
    return Av/101

def Diff(jet,label, y, k, i ):
    return (si(jet, label, y, k, i) - Av(jet, label, y, k))**2

def STDdev(jet, label, y, k):
    s = np.zeros((y), dtype=np.float64)
    for i in range(101):
        s = Diff(jet, label, y, k, i) + s
    return np.sqrt(s/100)

def Cpdf(jet,label,y,k):
    Cpdf = np.zeros((y,y), dtype = np.float64)
    for i in range(y):
        for j in range(y):
            Cpdf[i,j] = STDdev(jet,label,y,k)[i] * STDdev(jet,label,y,k)[j] #outer product of the errors to find the Cpdf
    return Cpdf

def shelve(jet, label, y):

    for k in range(121,122): #(114-125 for NN23, 115, 117, 118, 119, 121 for NN30
        values = Cpdf(jet,label,y,k)
        Cpd4['NN30_0'+ str(k)+''] = values 
    return Cpd4  
    
if __name__ =="__main__":
    y2 = 13
    y3 = 13
    y4 = 12  

    shelve(theory_z4, 'NNPDF30_0', y4)

theory_z2.close()
theory_z3.close()
theory_z4.close()
Cpd2.close()
Cpd3.close()
Cpd4.close()