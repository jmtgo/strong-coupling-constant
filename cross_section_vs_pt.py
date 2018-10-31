'''this module plots a graph of the measured cross section as a function of the leading jet
for inclusive z+ >_1,2,3,4 jet events, and the error as a fraction of one.'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shelve
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
 
data_z1 = pd.read_csv('HEPData-ins1514251-v2-Table_15.csv') #loading the data
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
theory_z2 = shelve.open('2jet')
theory_z3 = shelve.open('threejet')
theory_z4 = shelve.open('fourjet')

font = {'family': 'serif', #setting the font size
        'weight': 'normal',
        'size': 12,}
font1 = {'family': 'serif',
        'weight': 'normal',
        'size': 10,}


def CHI_SQ(st, l, sy, y, yth):
    C = np.diag(quad(st,l,sy)**2)
    X = np.zeros(16, dtype = np.float64)   
    for i in np.arange(112, 128): 
        X[i-112] = np.inner((y - np.array(yth['0'+str(i)+''])[1:len(y)+1] ), (np.inner(np.linalg.inv(C),(y - np.array(yth['0'+str(i)+''])[1:len(y)+1] ).T)))
    return X

def CHI_graph(ax, st2, l2, sy2, y2, yth2, st3, l3, sy3, y3, yth3, st4, l4, sy4, y4, yth4):
    fig1 =plt.figure(1)
   
    plt.scatter(ax,CHI_SQ(st2, l2, sy2, y2, yth2), marker = 'o', color = 'k')
    plt.scatter(ax,CHI_SQ(st3, l3, sy3, y3, yth3), marker = 'o', color = 'k')
    plt.scatter(ax,(CHI_SQ(st4, l4, sy4, y4, yth4))-350, marker = 'o', color = 'k')
    plt.plot(ax,CHI_SQ(st2, l2, sy2, y2, yth2), label = r'${\rm 2 \ jets}$' , color = 'b')
    plt.plot(ax,CHI_SQ(st3, l3, sy3, y3, yth3), label = r'${\rm 3 \ jets}$', color='r')
    plt.plot(ax,CHI_SQ(st4, l4, sy4, y4, yth4)-350, label = r'${\rm 4 \ jets,\ \chi ^2 \ -350}$', color='g')
    plt.xlabel(r'${\rm \alpha _s}$')
    plt.ylabel(r'${\chi ^2(\alpha_s)}$')
    plt.legend(loc = 'upper left')
    plt.show()

def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='k', edgecolor='None', alpha=0.25):

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
    ax1 = fig1.add_axes((.1,.44,.8,.52,), label = 'axes1')
    ax1.scatter(x1,y1, c='k', marker = 'o')
    ax1.scatter(x2,10e-2*y2, c='k', marker = 'o')
    ax1.scatter(x3,10e-3*y3, c='k', marker = 'o')
    ax1.scatter(x4,10e-4*y4, c='k', marker = 'o')
    plt.errorbar(x1,y1,yerr=None, xerr = data_z1['xerr'], fmt=u'None') #plotting x error bars
    plt.errorbar(x2,10e-2*y2,yerr=None, xerr = data_z2['xerr'], fmt=u'None')
    plt.errorbar(x3,10e-3*y3,yerr=None, xerr = data_z3['xerr'], fmt=u'None')
    plt.errorbar(x4,10e-4*y4,yerr=None, xerr = data_z4['xerr'], fmt=u'None')
    ax1.set_yscale('log') #makes the y scale logarithmic
    plt.ylim(20e-9,10e1)
    plt.xlim(20,700)
    ax1.get_xaxis().set_ticklabels([]) #no x tick labels wanted
    ax1.xaxis.set_minor_locator(plt.MaxNLocator(34)) #sets how many markers on the x axis
    plt.ylabel(r'${\rm d\sigma /dp_{\rm T}^{\rm jet} \, [pb/GeV]}$', fontdict= font)
    plt.text(100, 10, r'${\rm ATLAS\, 13\, TeV,\, 3.16 fb ^{\rm -1}}$', fontdict=font) #adding labels to the curves
    plt.text(400,10, r'${\rm Z/ \gamma^\star ( \rightarrow l^+l^- ) + \geq 1 \, jet }$')
    plt.text(100,1, r'${\rm Z/ \gamma^\star + \geq 1 \, jet}$', fontdict=font1, rotation=-20)
    plt.text(100,4e-2, r'${\rm Z/ \gamma^\star + \geq 2\,jets \, x10^{\rm -1}}$', fontdict=font1, rotation=-15)
    plt.text(100,10e-4, r'${\rm Z/ \gamma^\star + \geq 3\,jets \, x10^{\rm -2}}$', fontdict=font1, rotation=-10)
    plt.text(100,4e-5,r'${\rm Z/ \gamma^\star + \geq 4\,jets \, x10^{\rm -3}}$', fontdict=font1, rotation=-10)
    plt.text(40, 10e-8, r'${\rm anti-k_t \, jets, R=0.4,\, p_{\rm T}^{\rm jet} \rm > 30GeV, \,|y^{\rmjet}|<2.5}$', fontdict = font)

    ax2 = fig1.add_axes((.1,.35,.8,.09)) #adding error bar plots underneath the graph
    ax2.text(200,1.2, r'${\rm Z/ \gamma^\star + \geq 1 \,\, jet}$', fontdict=font1)
    ax2.plot(np.array([20.0,700]), np.array([1.0,1.0]), '-k')
    ax2.scatter(x1, np.ones(len(x1)), marker ='None')
    ax2.errorbar(x1, np.ones(len(x1)), yerr=yerror(y1, st1)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0) #statistical error lines
    plt.ylim(0.3,1.7)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax2.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax2.get_xaxis().set_ticklabels([])
    ax2.xaxis.set_minor_locator(plt.MaxNLocator(34))
    ax2.plot(np.array([20.0,700]), np.array([1.0,1.0]), '-k')
    ax2.scatter(x1, np.ones(len(x1)), marker ='None')
    _ = make_error_boxes(ax2, x1, np.ones(len(x1)), xerror(np.array(data_z1['xerr'])), yerror(y1, quad(st1, l1, sy1))) #calling error box function to draw shaded boxes
    


    ax3 = fig1.add_axes((.1,.26,.8,.09))
    ax3.text(200,1.2, r'${\rm Z/ \gamma^\star + \geq 2\,\,jets}$', fontdict=font1)
    plt.ylabel(r'$\rm Pred./Data$')
    ax3.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax3.scatter(x2, np.ones(len(x2)), marker ='None')
    ax3.errorbar(x2, np.ones(len(x2)), yerr=yerror(y2, st2)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0) #statistcal error lines, plots total yerror
    plt.ylim(0.3,1.7)
    ax3.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax3.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax3.get_xaxis().set_ticklabels([])
    ax3.xaxis.set_minor_locator(plt.MaxNLocator(34))
    ax3.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax3.scatter(x2, np.ones(len(x2)), marker ='None')
    ax3.scatter(x2, (np.array(theory_z2['0118'])[1:len(np.array(theory_z2['0118']))-1]/y2), marker = 'x' )
    _ = make_error_boxes(ax3, x2, np.ones(len(x2)), xerror(np.array(data_z2['xerr'])), yerror(y2, quad(st2, l2, sy2)))

    ax4 = fig1.add_axes((.1,.17,.8,.09))
    ax4.text(200,1.2, r'${\rm Z/ \gamma^\star + \geq 3\,\,jets}$', fontdict=font1)
    ax4.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax4.scatter(x3, np.ones(len(x3)), marker ='None')
    ax4.errorbar(x3, np.ones(len(x3)), yerr=yerror(y3, st3)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0) #statistical error lines
    plt.ylim(0.2,1.7)
    ax4.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax4.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax4.get_xaxis().set_ticklabels([])
    ax4.xaxis.set_minor_locator(plt.MaxNLocator(34))
    ax4.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax4.scatter(x3, np.ones(len(x3)), marker ='None')
    ax4.scatter(x3, (np.array(theory_z3['0118'])[1:len(np.array(theory_z3['0118']))-1]/y3), marker = 'x' )
    _ = make_error_boxes(ax4, x3, np.ones(len(x3)), xerror(np.array(data_z3['xerr'])), yerror(y3, quad(st3, l3, sy3)))



    ax5 = fig1.add_axes((.1,.09, .8,.09))
    ax5.text(200,1.2, r'${\rm Z/ \gamma^\star + \geq 4\,\,jets}$', fontdict=font1)
    ax5.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax5.scatter(x4, np.ones(len(x4)), marker ='None')
    ax5.errorbar(x4, np.ones(len(x4)), yerr=yerror(y4, st4)[0,:], xerr =None, fmt=u'None', ecolor= 'b', capsize=0) #statistical error lines
    plt.ylim(0.3,1.7)
    ax5.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax5.yaxis.set_minor_locator(plt.MaxNLocator(4))
    plt.yticks(np.arange(0.5,2.0,0.5, dtype = float), color = 'k')
    plt.xlim(20,700)
    ax5.xaxis.set_minor_locator(plt.MaxNLocator(34))
    ax5.plot(np.array([20.0,500]), np.array([1.0,1.0]), '-k')
    ax5.scatter(x4, np.ones(len(x4)), marker ='None')
    ax5.scatter(x4, (np.array(theory_z4['0118'])[1:len(np.array(theory_z4['0118']))-2]/y4), marker = 'x' )
    _ = make_error_boxes(ax5, x4, np.ones(len(x4)), xerror(np.array(data_z4['xerr'])), yerror(y4, quad(st4, l4, sy4)))
    plt.xlabel(r'$\rm p_{\rm T}^{\rm jet} \,(leading\, jet) \,[GeV]}$', fontdict = font)
    return plt.show()


def yerror(y,e): #how to calculate the errors to be plotted as an (n,2) matrix for error box
    yerr = np.zeros((2, len(y)))
    yerr[0,:] = ((y+e)/y)-1
    yerr[1,:] = ((y+e)/y)-1
    return yerr

def xerror(xe): #making the xerrors into (n,2) matrix for error box
    xerror = np.zeros((2, len(xe)))
    xerror[0,:] = xe
    xerror[1,:] = xe
    return xerror

def quad(st, l, sy): #calulating total errors in y
    return np.sqrt(st**2 +l**2 +sy**2)

if __name__ =="__main__": 
    x1 = np.array(data_z1['pT(jet) [GeV]'])
    x2 = np.array(data_z2['pT(jet) [GeV]'])
    x3 = np.array(data_z3['pT(jet) [GeV]'])
    x4 = np.array(data_z4['pT(jet) [GeV]'])
    y1 = np.array(data_z1['Cross_Section_Zll [pb]'])
    y2 = np.array(data_z2['Cross_Section_Zll [pb]'])
    y3 = np.array(data_z3['Cross_Section_Zll [pb]'])
    y4 = np.array(data_z4['Cross_Section_Zll [pb]'])
    st1 = np.array(data_z1['stat +']) #statistical error
    st2 = np.array(data_z2['stat +'])
    st3 = np.array(data_z3['stat +'])
    st4 = np.array(data_z4['stat +'])
    sy1 = np.array(data_z1['syst +']) #systematic error
    sy2 = np.array(data_z2['syst +'])
    sy3 = np.array(data_z3['syst +'])
    sy4 = np.array(data_z4['syst +'])
    l1 = np.array(data_z1['lumi +']) #lumi error
    l2 = np.array(data_z2['lumi +'])
    l3 = np.array(data_z3['lumi +'])
    l4 = np.array(data_z4['lumi +'])
    ax = np.arange(0.112,0.128, 0.001)

    graph = graph(x1,x2,x3,x4,y1,y2,y3,y4, st1, st2,st3, st4, l1, l2, l3, l4, sy1, sy2, sy3, sy4)
    CHI_graph(ax, st2, l2, sy2, y2, theory_z2, st3, l3, sy3, y3, theory_z3, st4, l4, sy4, y4, theory_z4)
  
theory_z2.close()
theory_z3.close()
theory_z4.close()
