'''this module plots a graph of the measured cross section as a function of the leading jet
for inclusive z+ >_1,2,3,4 jet events'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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


def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='k',
                     edgecolor='None', alpha=0.25):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata.T, ydata.T, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha, 
                         edgecolor=edgecolor)
    # Add collection to axes
    artists = ax.add_collection(pc)
    # Plot errorbars
    #artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,fmt='None', ecolor='None')
    return artists

def graph(x1,x2,x3,x4,y1,y2,y3,y4, st1, st2,st3, st4, l1, l2, l3, l4, sy1, sy2, sy3, sy4):
    fig1 =plt.figure(1)
    #ax1 =fig.add_subplot(211)
    ax1 = fig1.add_axes((.1,.36,.8,.6,), label = 'axes1')
    ax1.scatter(x1,y1, c='k', marker = 'o')
    ax1.scatter(x2,y2, c='k', marker = 'o')
    ax1.scatter(x3,y3, c='k', marker = 'o')
    ax1.scatter(x4,y4, c='k', marker = 'o')
    plt.errorbar(x1,y1,yerr=None, xerr = data_z1['xerr'], fmt=u'None')
    plt.errorbar(x2,y2,yerr=None, xerr = data_z2['xerr'], fmt=u'None')
    plt.errorbar(x3,y3,yerr=None, xerr = data_z3['xerr'], fmt=u'None')
    plt.errorbar(x4,y4,yerr=None, xerr = data_z4['xerr'], fmt=u'None')
    ax1.set_yscale('log')
    plt.ylim(20e-9,10e1)
    plt.xlim(20,700)
    ax1.get_xaxis().set_ticklabels([])
    ax1.xaxis.set_minor_locator(plt.MaxNLocator(34))
    plt.ylabel('d'r'${\rm \sigma}$''\d'r'$\rm p_{\rm T}^{\rm jet}$' ' [pb/GeV]', fontdict= font)
    plt.text(100, 10, 'ATLAS 13 TeV, 3.16 fb' r'$^{\rm -1}$', fontdict=font)
    plt.text(400,10, r'${\rm Z/ \gamma^\star ( \rightarrow l^+l^- ) + \geq 1 jet }$')
    plt.text(100,1, r'${\rm Z/ \gamma^\star + \geq 1 jet}$', fontdict=font1, rotation=-30)
    plt.text(100,4e-2, r'${\rm Z/ \gamma^\star + \geq}$' '2 jets,' r'$\rm x10^{\rm -1}$', fontdict=font1, rotation=-25)
    plt.text(100,10e-4, r'${\rm Z/ \gamma^\star + \geq}$' '3 jets, ' r'$\rm x10^{\rm -2}$', fontdict=font1, rotation=-20)
    plt.text(100,4e-5, r'${\rm Z/ \gamma^\star + \geq}$' '4 jets, ' r'$\rm x10^{\rm -3}$', fontdict=font1, rotation=-20)

    ax2 = fig1.add_axes((.1,.27,.8,.09))
    ax2.plot(np.array([20.0,700]), np.array([1.0,1.0]), '-k')
    ax2.scatter(x1, np.ones(len(x1)), marker ='None')
    ax2.errorbar(x1, np.ones(len(x1)), yerr=yerror(y1, st1)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0)
    plt.ylim(0.5,1.5)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax2.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax2.get_xaxis().set_ticklabels([])
    ax2.xaxis.set_minor_locator(plt.MaxNLocator(34))
    ax2.plot(np.array([20.0,700]), np.array([1.0,1.0]), '-k')
    ax2.scatter(x1, np.ones(len(x1)), marker ='None')
    _ = make_error_boxes(ax2, x1, np.ones(len(x1)), xerror(np.array(data_z1['xerr'])), yerror(y1, quad(st1, l1, sy1)))
    


    ax3 = fig1.add_axes((.1,.18,.8,.09))
    ax3.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax3.scatter(x2, np.ones(len(x2)), marker ='None')
    ax3.errorbar(x2, np.ones(len(x2)), yerr=yerror(y2, st2)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0)
    plt.ylim(0.5,1.5)
    ax3.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax3.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax3.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax3.scatter(x2, np.ones(len(x2)), marker ='None')
    _ = make_error_boxes(ax3, x2, np.ones(len(x2)), xerror(np.array(data_z2['xerr'])), yerror(y2, quad(st2, l2, sy2)))
    plt.xlabel(r'$\rm p_{\rm T}^{\rm jet}$'' (leading jet) [GeV]', fontdict = font)

    #ax4 = fig1.add_axes((.1,.09,.8,.09))
    #ax5 = fig1.add_axes((.1,.01,.8,.09))
    return plt.show()


def yerror(y,e):
    yerr = np.zeros((2, len(y)))
    yerr[0,:] = ((y+e)/y)-1
    yerr[1,:] = ((y+e)/y)-1
    return yerr

def xerror(xe):
    xerror = np.zeros((2, len(xe)))
    xerror[0,:] = xe
    xerror[1,:] = xe
    return xerror

def quad(st, l, sy):
    return np.sqrt(st**2 +l**2 +sy**2)

if __name__ =="__main__": 
    x1 = np.array(data_z1['pT(jet) [GeV]'])
    x2 = np.array(data_z2['pT(jet) [GeV]'])
    x3 = np.array(data_z3['pT(jet) [GeV]'])
    x4 = np.array(data_z4['pT(jet) [GeV]'])
    y1 = np.array(data_z1['Cross_Section_Zll [pb]'])
    y2 = np.array(10e-2*data_z2['Cross_Section_Zll [pb]'])
    y3 = np.array(10e-3*data_z3['Cross_Section_Zll [pb]'])
    y4 = np.array(10e-4* data_z4['Cross_Section_Zll [pb]'])
    st1 = np.array(data_z1['stat +'])
    st2 = np.array(data_z2['stat +'])
    st3 = np.array(data_z3['stat +'])
    st4 = np.array(data_z4['stat +'])
    sy1 = np.array(data_z1['syst +'])
    sy2 = np.array(data_z2['syst +'])
    sy3 = np.array(data_z3['syst +'])
    sy4 = np.array(data_z4['syst +'])
    l1 = np.array(data_z1['lumi +'])
    l2 = np.array(data_z2['lumi +'])
    l3 = np.array(data_z3['lumi +'])
    l4 = np.array(data_z4['lumi +'])
    graph = graph(x1,x2,x3,x4,y1,y2,y3,y4, st1, st2,st3, st4, l1, l2, l3, l4, sy1, sy2, sy3, sy4)
