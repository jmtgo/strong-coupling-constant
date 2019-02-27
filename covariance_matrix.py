'''This module is designed to evaluate the covariance error matrix for multiple jets.'''
import numpy as np
import pandas as pd
import shelve
import matplotlib.pyplot as plt

Cov_mega = shelve.open('Cov_all')
#load the error excel sheets for electron channel
error_2es = pd.read_csv('2jetsysterrors_elec.csv') 
error_3es = pd.read_csv('HEPData-ins1514251-v2-Table_49.csv')
error_4es = pd.read_csv('HEPData-ins1514251-v2-Table_51.csv')

#load the error excel sheets for muon channel
error_2us = pd.read_csv('HEPData-ins1514251-v2-Table_48.csv')
error_3us = pd.read_csv('HEPData-ins1514251-v2-Table_50.csv')
error_4us = pd.read_csv('HEPData-ins1514251-v2-Table_52.csv')


#load the cross section for each jet, for each channel
exp_2e = pd.read_csv('HEPData-ins1514251-v2-Table_16.csv') #electron
exp_3e = pd.read_csv('HEPData-ins1514251-v2-Table_19.csv')
exp_4e = pd.read_csv('HEPData-ins1514251-v2-Table_22.csv')

exp_2u = pd.read_csv('HEPData-ins1514251-v2-Table_17.csv') #muon
exp_3u = pd.read_csv('HEPData-ins1514251-v2-Table_20.csv')
exp_4u = pd.read_csv('HEPData-ins1514251-v2-Table_23.csv')


def full_error_vec(er_2e,er_3e, er_4e, er_2u, er_3u, er_4u): #makes the giant matirx with all the systematic and luminosity errors on.
    A = np.concatenate((er_2e,er_3e), axis=1)
    B = np.concatenate((A,er_4e), axis=1)
    C = np.concatenate((B,er_2u), axis=1)
    D = np.concatenate((C,er_3u), axis=1)
    E = np.concatenate((D,er_4u), axis=1)
    return E


def full_cross_vec(xsec_2e,xsec_3e,xsec_4e,xsec_2u,xsec_3u,xsec_4u): #makes the full vector with all the channels and jets in.
    A = np.append(xsec_2e,xsec_3e)
    B = np.append(A, xsec_4e)
    C = np.append(B, xsec_2u)
    D = np.append(C, xsec_3u)
    E = np.append(D, xsec_4u)
    print np.shape(E)
    return E

def full_stat_vec(st_e2,st_e3,st_e4,st_u2,st_u3,st_u4):
    A = np.append(st_e2,st_e3)
    B = np.append(A, st_e4)
    C = np.append(B, st_u2)
    D = np.append(C, st_u3)
    E = np.append(D, st_u4)
    return E

def Cov_s(e): 
    # for each source of systematic error a covariance matrix is calculated via e[i] * e[j]
    A = np.zeros((len(e),len(e)), dtype = np.float64)
    for i in np.arange(len(e)):
        for j in np.arange(len(e)):
            A[i,j]= e[i]*e[j]
    return A

def Cov_syst2(e): 
    # all the systematic error types matrices are then added together
    Csy2 = np.zeros((len(e[0,:]),len(e[0,:])), dtype = np.float64)
    for k in np.arange(len(e[:,0])):
        Csy2 = Cov_s(e[k,:]) + Csy2
    return Csy2

def Cov_syst(e):
    #takes the square-root of all the elements in the covariance matrix as a percentage of the measurement
    return np.sqrt(Cov_syst2(e))

def Cov_syst100(e, y):
    #finds what the actual error is by "reverse percentaging" of the actual measurment
    X = np.zeros((len(e[0,:]),len(e[0,:])), dtype= np.float64)
    for i in np.arange(len(e[0,:])):
        for j in np.arange(len(e[0,:])):
            print i,j
            X[i,j] = ((Cov_syst(e)[i,j])/100 * np.sqrt(y[i]*y[j]))**2
    return X

def Cov_sylu(e,y):
    #squares all the elements in the covariance matrix as they should have those dimensions
    X = np.zeros((len(e[0,:]),len(e[0,:])), dtype= np.float64)
    for i in np.arange(len(e[0,:])):
        for j in np.arange(len(e[0,:])):
            print 'second set', i,j
            X[i,j] = Cov_syst100(e,y)[i,j] **2
    return X


def Cov_stat(st):
    #calculates the statistcal error matrix 
    return np.diag(st**2)

def Cov_full(st, Cov_SYLU):
    #calculates the full covariance matrix
    return Cov_SYLU + Cov_stat(st)

if __name__ =="__main__": 
    #electron error tables
    er_2e = np.array(error_2es.iloc[1:,1:])
    er_3e = np.array(error_3es.iloc[1:,1:])
    er_4e = np.array(error_4es.iloc[1:,1:])

    #muon error tables
    er_2u = np.array(error_2us.iloc[1:,1:])
    er_3u = np.array(error_3us.iloc[1:,1:])
    er_4u = np.array(error_4us.iloc[1:,1:])
    
    #cross section data for electron channel
    xsec_2e = np.array(exp_2e['Cross_Section_Zee [pb]'])
    xsec_3e = np.array(exp_3e['Cross_Section_Zee [pb]'])
    xsec_4e = np.array(exp_4e['Cross_Section_Zee [pb]'])
    
    #cross section data for muon channel
    xsec_2u = np.array(exp_2u['Cross_Section_Zmm [pb]'])
    xsec_3u = np.array(exp_3u['Cross_Section_Zmm [pb]'])
    xsec_4u = np.array(exp_4u['Cross_Section_Zmm [pb]'])

    #statistical errors for electron channel
    st_e2 = np.array(exp_2e['stat +'])
    st_e3 = np.array(exp_3e['stat +'])
    st_e4 = np.array(exp_4e['stat +'])

    #statistical errors for muon channel
    st_u2 = np.array(exp_2u['stat +'])
    st_u3 = np.array(exp_3u['stat +'])
    st_u4 = np.array(exp_4u['stat +'])

    #systematic errors for electron channel
    sy_e2 = np.array(exp_2e['syst +'])
    sy_e3 = np.array(exp_3e['syst +'])
    sy_e4 = np.array(exp_4e['syst +'])

    #systematic errors for muon channel
    sy_u2 = np.array(exp_2u['syst +'])
    sy_u3 = np.array(exp_3u['syst +'])
    sy_u4 = np.array(exp_4u['syst +'])

    #luminosity errors for electron channel
    lu_e2 = np.array(exp_2e['lumi +'])
    lu_e3 = np.array(exp_3e['lumi +'])
    lu_e4 = np.array(exp_4e['lumi +'])

    #luminosity errors for muon channel
    lu_u2 = np.array(exp_2u['lumi +'])
    lu_u3 = np.array(exp_3u['lumi +'])
    lu_u4 = np.array(exp_4u['lumi +'])


    


    syst_Errors = full_stat_vec(sy_e2,sy_e3,sy_e4,sy_u2,sy_u3,sy_u4)
    lumi_Errors = full_stat_vec(lu_e2,lu_e3,lu_e4,lu_u2,lu_u3,lu_u4)
    Total = syst_Errors**2 + lumi_Errors**2
    #print np.shape(er_2e)
    #print np.shape(er_4e)
    Cov_mega['FULL_error_vec'] = full_error_vec(er_2e,er_3e, er_4e, er_2u, er_3u, er_4u)
    Cov_mega['FULL_xsec_vec'] = full_cross_vec(xsec_2e,xsec_3e,xsec_4e,xsec_2u,xsec_3u,xsec_4u)
    #full_cross_vec(xsec_2e,xsec_3e,xsec_4e,xsec_2u,xsec_3u,xsec_4u)
    #Cov_mega['FULL_SYST+LUMI'] = Cov_syst100(full_error_vec(er_2e,er_3e, er_4e, er_2u, er_3u, er_4u),full_cross_vec(xsec_2e,xsec_3e,xsec_4e,xsec_2u,xsec_3u,xsec_4u))
    #Cov_mega['FULL_ALLTHREE'] = Cov_full(full_stat_vec(st_e2,st_e3,st_e4,st_u2,st_u3, st_u4),Cov_mega['FULL_SYST+LUMI'])
    #print np.diag(Cov_mega['FULL_ALLTHREE'])
    print np.max(np.sqrt(np.diag(Cov_mega['FULL_SYST+LUMI']))/ np.sqrt(Total))
    print np.min(np.sqrt(np.diag(Cov_mega['FULL_SYST+LUMI']))/ np.sqrt(Total))
Cov_mega.close()


    


