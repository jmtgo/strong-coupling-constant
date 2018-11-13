'''This module is designed to evaluate the covariance error matrix for multiple jets, initially for the systematic and lumi errors only'''

import numpy as np
import pandas as pd
import shelve

raw_2e = pd.read_csv('HEPData-ins1514251-v2-Table_47_2.csv') #_2 indicated that the muon errors are now included 
raw_3e = pd.read_csv('HEPData-ins1514251-v2-Table_49.csv')
raw_4e = pd.read_csv('HEPData-ins1514251-v2-Table_51.csv')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
Cov_all = shelve.open('Cov_all')

def Cov_s(e): 
    # for each source of systematic error a covariance matrix is calculated via e[i] * e[j]
    A = np.zeros((len(e),len(e)), dtype = np.float64)
    for i in np.arange(len(e)):
        for j in np.arange(len(e)):
            A[i,j]= e[i]*e[j]
    return A

def Cov_syst2(e): 
    # all the systematic error matrices are then added together
    Csy2 = np.zeros((len(e[0,:]),len(e[0,:])), dtype = np.float64)
    for k in np.arange(len(e[:,0])):
        Csy2 = Cov_s(e[k,:]) + Csy2
    return Csy2

def Cov_syst(e):
    #takes the square-root of all the elements in the covariance matrix as a percentage of the measurement
    return np.sqrt(Cov_syst2(e))

def Cov_syst100(e, y):
    #finds what the error is by "reverse percentaging" of the actual measurment
    X = np.zeros((len(e[0,:]),len(e[0,:])), dtype= np.float64)
    for i in np.arange(len(e[0,:])):
        for j in np.arange(len(e[0,:])):
            X[i,j] = (Cov_syst(e)[i,j])/100 * np.sqrt(y[i]*y[j])
    return X

def Cov_sylu(e,y):
    #squares all the elements in the covariance matrix as they should have those dimensions
    X = np.zeros((len(e[0,:]),len(e[0,:])), dtype= np.float64)
    for i in np.arange(len(e[0,:])):
        for j in np.arange(len(e[0,:])):
            X[i,j] = Cov_syst100(e,y)[i,j] **2
    return X


def Cov_stat(st):
    #calculates the statistcal error matrix 
    return np.diag(st**2)

def Cov_full(e, st,y):
    #calculates the full covariance matrix
    return Cov_sylu(e,y) + Cov_stat(st)

if __name__ =="__main__": 

    twoe = np.array(raw_2e.iloc[0:,1:])
    st2 = np.array(data_z2['stat +'])
    three = np.array(raw_3e.iloc[0:,1:])
    st3 = np.array(data_z3['stat +'])
    foure = np.array(raw_4e.iloc[0:,1:])
    st4 = np.array(data_z4['stat +'])
    y2 = np.array(data_z2['Cross_Section_Zll [pb]'])
    y3 = np.array(data_z3['Cross_Section_Zll [pb]'])
    y4 = np.array(data_z4['Cross_Section_Zll [pb]'])


    Cov_all['twoj'] = Cov_full(twoe,st2, y2)
    Cov_all['threej'] = Cov_full(three, st3, y3)
    Cov_all['fourj'] = Cov_full(foure,st4, y4)

Cov_all.close()

    


