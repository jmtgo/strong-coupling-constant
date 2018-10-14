'''this module plots a graph of the measured cross section as a function of the leading jet
for inclusive z+ >_1,2,3,4 jet events'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data_z1 = pd.read_csv('HEPData-ins1514251-v2-Table_15.csv')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')

font = {'family': 'serif',
        'weight': 'normal',
        'size': 12,}
font1 = {'family': 'serif',
        'weight': 'normal',
        'size': 10,}

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
    plt.xlabel(r'$\rm p_{\rm T}^{\rm jet}$'' (leading jet) [GeV]', fontdict = font)
    plt.ylabel('d'r'${\rm \sigma}$''\d'r'$\rm p_{\rm T}^{\rm jet}$' ' [pb/GeV]', fontdict= font)
    plt.text(100, 10, 'ATLAS 13 TeV, 3.16 fb' r'$^{\rm -1}$', fontdict=font)
    plt.text(400,10, r'${\rm Z/ \gamma^\star (->l^+l^- ) + \geq 1 jet }$')
    plt.text(100,1, r'${\rm Z/ \gamma^\star + \geq 1 jet}$', fontdict=font1, rotation=-30)
    plt.text(100,4e-2, r'${\rm Z/ \gamma^\star + \geq}$' '2 jets,' r'$\rm x10^{\rm -1}$', fontdict=font1, rotation=-25)
    plt.text(100,10e-4, r'${\rm Z/ \gamma^\star + \geq}$' '3 jets, ' r'$\rm x10^{\rm -2}$', fontdict=font1, rotation=-20)
    plt.text(100,4e-5, r'${\rm Z/ \gamma^\star + \geq}$' '4 jets, ' r'$\rm x10^{\rm -3}$', fontdict=font1, rotation=-20)
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