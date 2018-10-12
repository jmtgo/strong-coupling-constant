'''this module plots a graph of the measured cross section as a function of the leading jet
for inclusive z+ >_1,2,3,4 jet events'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data_z1 = pd.read_csv('HEPData-ins1514251-v2-Table_15.csv')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')

#print np.shape(data_z4['Cross_Section_Zll [pb]'])
#print (np.shape(data_z4['pT(jet) [GeV]']))

def graph(x1,x2,x3,x4,y1,y2,y3,y4):
    fig =plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x1,y1, c='k', marker = 'o')
    ax1.scatter(x2,y2, c='k', marker = 'o')
    ax1.scatter(x3,y3, c='k', marker = 'o')
    ax1.scatter(x4,y4, c='k', marker = 'o')
    plt.errorbar(x1,y1,yerr=None, xerr = data_z1['xerr'], fmt=u'None')
    plt.errorbar(x2,y2,yerr=None, xerr = data_z2['xerr'], fmt=u'None')
    plt.errorbar(x3,y3,yerr=None, xerr = data_z3['xerr'], fmt=u'None')
    plt.errorbar(x4,y4,yerr=None, xerr = data_z4['xerr'], fmt=u'None')
    ax1.set_yscale('log')
    plt.ylim(10e-8,10e1)
    plt.xlim(20,700)
    plt.xlabel('p,,T,,^jet^')
    return plt.show()

if __name__ =="__main__": 
    x1 = data_z1['pT(jet) [GeV]']
    x2 = data_z2['pT(jet) [GeV]']
    x3 = data_z3['pT(jet) [GeV]']
    x4 = data_z4['pT(jet) [GeV]']
    y1 = data_z1['Cross_Section_Zll [pb]']
    y2 = 10e-2*data_z2['Cross_Section_Zll [pb]']
    y3 = 10e-3*data_z3['Cross_Section_Zll [pb]']
    y4 = 10e-4* data_z4['Cross_Section_Zll [pb]']
    graph = graph(x1,x2,x3,x4,y1,y2,y3,y4)
    plt.show()