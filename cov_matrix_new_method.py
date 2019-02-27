'''This module is designed to evaluate the covariance error matrix for multiple jets.'''
import numpy as np
import pandas as pd
import shelve
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal

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


def full_stat_vec(st_e2,st_e3,st_u2,st_u3):
    A = np.append(st_e2,st_e3)
    C = np.append(A, st_u2)
    D = np.append(C, st_u3)
    return D

def re_scale_error(e,y):
    E = np.zeros((len(e[:,0]), len(e[0,:])))
    for i in range(len(e[0,:])):
        E[:,i] = y[i] * e[:,i]/100
    return E
    
def Cov_s(e,y,k): 
    # for each source of systematic error a covariance matrix is calculated via e[i] * e[j]
    A = np.zeros((len(e[0,:]),len(e[0,:])), dtype = np.float64)
    for i in np.arange(len(e[0,:])):
        for j in np.arange(len(e[0,:])):
            A[i,j]= re_scale_error(e,y)[k,i]*re_scale_error(e,y)[k,j]
    return A

def Cov_syst2(e,y): 
    # all the systematic error types matrices are then added together
    Csy2 = np.zeros((len(e[0,:]),len(e[0,:])), dtype = np.float64)
    for k in np.arange(len(e[:,0])):
        print k
        Csy2 = Cov_s(e,y,k) + Csy2
    return Csy2

def total_error_other(e, bin,y):
    return np.matmul(bin,np.matmul(Cov_syst2(e,y)[13:26,13:26],bin.T)) 

def Cov_stat(st):
    #calculates the statistcal error matrix 
    return np.diag(st**2)

def Cov_full(st, Cov_SYLU):
    #calculates the full covariance matrix
    return Cov_SYLU + Cov_stat(st)

def colourplot(x,y, array):

    plt.subplot(111)
    pcm = plt.pcolor(x,y,array, vmin = np.amin(array) ,vmax= np.amax(array),cmap= 'Purples')
    cbar = plt.colorbar(pcm)
    cbar.set_label('value of element of Cov_matrix')
    plt.show()

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

    #bin vector
    binx= np.array([20,20,20,20,20,20,20,30,50,50,50,50,100])
    binx4= np.array([20,20,20,20,20,20,20,30,50,50,100,100])
    binx_full = np.append(binx,np.append(binx, np.append(binx,binx)))
    e = Cov_mega['FULL_error_vec_no4']
    y = Cov_mega['FULL_xsec_vec_no4']
    

    #print np.count_nonzero(re_scale_error(Cov_mega['FULL_error_vec_no4'], Cov_mega['FULL_xsec_vec_no4'] ))
    #print np.shape(re_scale_error(Cov_mega['FULL_error_vec_no4'], Cov_mega['FULL_xsec_vec_no4'] )[0,:])


    #print np.count_nonzero(Cov_syst2(e,y)[0:12,0:12])
    yc = np.arange(52)
    x = np.arange(52)

    array = Cov_syst2(e,y)

    print np.amin(array)
    colourplot(x,yc, array)

    



    #Cov_mega['FULL_error_vec_no4'] = full_error_vec(er_2e,er_3e, er_2u, er_3u,)
    #Cov_mega['FULL_xsec_vec_no4'] = full_cross_vec(xsec_2e,xsec_3e,xsec_2u,xsec_3u)
Cov_mega.close()


    


