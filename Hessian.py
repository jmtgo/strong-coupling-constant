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

def si(jet, label,y,k):
    return np.array(jet[label + str(k) + ''])[1:y+1]

def Up_Diff(jet,label, y, k):
    D = np.amax([np.zeros(y), (si(jet, label, y, k+1)- si(jet,label, y, 0)), (si(jet, label, y, k)- si(jet,label, y, 0))], axis=0)
    return D

def UpE_sq(jet, label, y):
    UE = np.zeros(y)
    for k in range(1,53,2): #53 for CT10, 57 for CT14
        UE = (Up_Diff(jet, label, y, k))**2 + UE
    return UE 

def Lo_Diff(jet,label, y, k):
    D = np.amin([np.zeros(y), (si(jet, label, y, k+1)+ si(jet,label, y, 0)), (si(jet, label, y, k)+ si(jet,label, y, 0))], axis=0)
    return D

def LoE_sq(jet, label, y):
    LE = np.zeros(y)
    for k in range(1,53,2): #53 for CT10, 57 for CT14
        LE = (Lo_Diff(jet, label, y, k))**2 + LE
    return LE 

def Av_e(jet, label, y):
    return (np.sqrt(LoE_sq(jet,label,y)) + np.sqrt(UpE_sq(jet, label, y)))/2

def Cpdf(jet,label,y):
    Cpdf = np.zeros((y,y), dtype = np.float64)
    for i in range(y):
        for j in range(y):
            Cpdf[i,j] = Av_e(jet,label,y)[i] * Av_e(jet,label,y)[j]
    return Cpdf

def Corr(jet,label,y):
    Corr = np.zeros((y,y), dtype=np.float64)
    for i in range(y):
        for j in range(y):
            Corr[i,j] = x
    return 

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
   
    print np.shape(Cpdf(theory_z2, 'CT10nlo', y2) )
   

theory_z2.close()
theory_z3.close()
theory_z4.close()