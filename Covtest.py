import shelve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
Cov = shelve.open('Cov_all')

def quad(st, l, sy): #calulating total errors in y
    return st**2 +l**2 +sy**2

def e_vs_a(label, st, l, sy ):
    C = np.array(Cov[label])
    R = np.zeros((len(st)), dtype = np.float64)
    for i in np.arange(len(st)):
        R[i] = C[i,i]/quad(st,l,sy)[i]
    return R

def R_graph(x,label, st, l, sy):
    plt.scatter (x, e_vs_a(label, st, l, sy))
    return plt.show


if __name__ =="__main__": 
   
    x2 = np.array(data_z2['pT(jet) [GeV]'])
    x3 = np.array(data_z3['pT(jet) [GeV]'])
    x4 = np.array(data_z4['pT(jet) [GeV]'])
 
    st2 = np.array(data_z2['stat +'])
    st3 = np.array(data_z3['stat +'])
    st4 = np.array(data_z4['stat +'])
   
    sy2 = np.array(data_z2['syst +'])
    sy3 = np.array(data_z3['syst +'])
    sy4 = np.array(data_z4['syst +'])

    l2 = np.array(data_z2['lumi +'])
    l3 = np.array(data_z3['lumi +'])
    l4 = np.array(data_z4['lumi +'])

    
    plt.show (R_graph(x2, 'twoj', st2, l2, sy2))
    plt.show(R_graph(x3, 'threej', st3, l3, sy3))
    plt.show(R_graph(x4, 'fourj', st4, l4, sy4))
Cov.close()