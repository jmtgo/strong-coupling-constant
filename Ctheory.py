'''This module evaluates the Ctheory matrix by slicing from 17x17 to the correct size by removing two zero rows and the higher pt bin rows'''

import shelve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


yth = shelve.open('theoryCM.shelve')
t2 = yth['Zee2j']
t3 = yth['Zee3j']
Ct = shelve.open('Ct')
Cp2 = shelve.open('Cpdf2')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
npc2_full = pd.read_csv('HEPData-ins1514251-v2-Table_66.csv')
npc3_full = pd.read_csv('HEPData-ins1514251-v2-Table_67.csv')


def si(jet,y): #gives the fastnlo cross section bin data.
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet)[2:y+2, 2:y+2] #for particular pdf parameter change k, is a vector length y

def np_matrix(npc):
    return np.diag(np.array(npc))

def C_Theory_matrix(npc, jet, y):
    return np.matmul(np_matrix(npc),np.matmul(si(jet,y), np_matrix(npc)))


if __name__ =="__main__": 

    y2 = 13
    y3 = 13
    y4 = 12
    y2e = np.array(data_z2['Cross_Section_Zll [pb]'])
    y3e = np.array(data_z3['Cross_Section_Zll [pb]'])
    y4e = np.array(data_z4['Cross_Section_Zll [pb]'])
    x= np.array([20,20,20,20,20,20,20,30,50,50,50,50,50,100])
    x4= np.array([20,20,20,20,20,20,20,30,50,50,50,100,100,])
    npc2 = npc2_full['NP_corrections_Zll [none]']
    npc3 = npc3_full['NP_corrections_Zll [none]']





    Ct['Ct3'] = C_Theory_matrix(npc3, t3, y3) 
    #print np.diag(t2)
    #print np.diag(si(t2,y2))
    #print np.diag(Ct['Ct2'])
    #print np.diag(si(t2,y2))

    #plotting the ratio of the diagonals of the covariance matrices for
    #CT14 Cpdf at principle a_s against Ctheory from DM.




    

yth.close()
Ct.close()
Cp2.close()