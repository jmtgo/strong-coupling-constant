'''This module evaluates the Cpdf covariance matrix for any pdf (not NN type)'''


import numpy as np
import shelve
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors

theory_z2 = shelve.open('twojet')
theory_z3 = shelve.open('threejet')
theory_z4 = shelve.open('fourjet')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
Cpd2 = shelve.open('Cpdf2')
Cpd3 = shelve.open('Cpdf3')
Cpd4 = shelve.open('Cpdf4')

def si(jet, label,y,k):
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet[label + str(k) + ''])[1:y+1] #for particular pdf parameter change k, is a vector length y

def Up_Diff(jet,label, y, k): #returns a vector length y.
    D = np.amax([np.zeros(y), (si(jet, label, y, k+1)- si(jet,label, y, 0)), (si(jet, label, y, k)- si(jet,label, y, 0))], axis=0)
    return D #finds the maximum of 0, si0-sik+, si0+sik-, where the k's come in pairs for plus and minus the parameter change. 

def UpE_sq(jet, label, y): #returns a vector
    UE = np.zeros(y)
    for k in range(1,41,2): #53 for CT10, 57 for CT14, 41 for MSTW, 51 for MMHT
        UE = (Up_Diff(jet, label, y, k))**2 + UE
    return UE #square all the D's for each k pair and adds them together.

def Lo_Diff(jet,label, y, k): #returns a vector length y
    D = np.amin([np.zeros(y), (si(jet, label, y, k+1)+ si(jet,label, y, 0)), (si(jet, label, y, k)+ si(jet,label, y, 0))], axis=0)
    return D #finds the minimum of 0, si0+sik+, si0 + sik-

def LoE_sq(jet, label, y): #returns a vector length y
    LE = np.zeros(y)
    for k in range(1,41,2): #53 for CT10, 57 for CT14, 41 for MSTW, 51 for MMHT
        LE = (Lo_Diff(jet, label, y, k))**2 + LE
    return LE #square all the D's for each k pair and adds them together.

def Av_e(jet, label, y): #averages the upper and lower errors. 
    return (np.sqrt(LoE_sq(jet,label,y)) + np.sqrt(UpE_sq(jet, label, y)))/2

def Cpdf(jet,label,y):
    Cpdf = np.zeros((y,y), dtype = np.float64)
    for i in range(y):
        for j in range(y):
            Cpdf[i,j] = Av_e(jet,label,y)[i] * Av_e(jet,label,y)[j] #outer product of the errors to find the Cpdf
    return Cpdf

def colormap(jet,label,y):
    X,Y = np.mgrid[0:y,0:y]
    pcm = plt.pcolor(X,Y, Cpdf(jet,label,y))
    return plt.show()



if __name__ =="__main__":
    x2 = np.array(data_z2['pT(jet) [GeV]'])
    x3 = np.array(data_z3['pT(jet) [GeV]'])
    x4 = np.array(data_z4['pT(jet) [GeV]'])
    y2 = 13
    y3 = 13
    y4 = 12   
   
    #Cpd2['MSTW08cl'] = Cpdf(theory_z2,'MSTW0868cl', y2)
    #Cpd3['MSTW08cl'] = Cpdf(theory_z3,'MSTW0868cl', y3)
    #Cpd4['MSTW08cl'] = Cpdf(theory_z4,'MSTW0868cl', y4)
  
   

theory_z2.close()
theory_z3.close()
theory_z4.close()
Cpd2.close()
Cpd3.close()
Cpd4.close()