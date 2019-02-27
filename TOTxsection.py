'''This module evaluates the total cross section for comparison of the predictions with the experimental data on a more overall level.
It also produces a ratio plot of the electtron channel systematic errors to theose given in HEP.'''

'''This module also evaluates and saves a figure showing the as dependance of the total cross section'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shelve
from scipy import optimize
from scipy.optimize import minimize

# these are the data for the ll channel with non broken down  systematic errors.

data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')

error_2es = pd.read_csv('2jetsysterrors_elec.csv') 
data_e2 = pd.read_csv('HEPData-ins1514251-v2-Table_16.csv')

Cov_all = shelve.open('Cov_all')


def tot_xsec_exp(y,binval):
    return np.sum(y* binval)

def scalevariation(binx, start, stop, step, pdf,n):
    #for a in [0.5,1,2]:
    #    for b in [0.5,1,2]:
    a = 1
    b = 1

    yt2 = shelve.open('pc2__('+str(a)+','+str(b)+')') #theoretical predictions
    yt3 = shelve.open('pc3__('+str(a)+','+str(b)+')')
    ytot = shelve.open('ytot')

    X = np.zeros((n), dtype= np.float64) #n is the number of alpha s values
    for i in range(start, stop, step):
        X[i] = tot_xsec_exp(yt3[pdf + str(i)],binx)
    ytot['three'] = X
    print ytot['three'] #mstw only 

    yt2.close()
    yt3.close()
    ytot.close()

def tot_Vs_as_plt(x):
    ytot = shelve.open('ytot')
    y2 = ytot['two'][1:22]
    y3 = y2 + ytot['three'][1:22] 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.scatter(x, y2, label='two jet inclusive cross section', color='k', marker = 'x')
    plt.scatter(x, y3, label = 'three jet inclusive cross section', color='k', marker = 'x')
    plt.plot(x, func(x, *curve_fit(x, y2)), label= 'two jet inclusive cross section')
    plt.plot(x, func(x, *curve_fit(x, y3)), label= 'three jet inclusive cross section')
    plt.xlim(0.104,0.13)
    plt.ylim(20,43)
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax.yaxis.set_minor_locator(plt.MaxNLocator(24))
    plt.yticks(np.arange(20,43,5, dtype = float), color = 'k')
    ax.vlines(x=0.12, ymin = 20, ymax = func(0.12, *curve_fit(x, y3)), colors='k', linestyles='solid')
    ax.vlines(x=0.115, ymin = 20, ymax = func(0.115, *curve_fit(x, y3)), colors='k', linestyles='solid')
    ax.hlines(y=func(0.12, *curve_fit(x, y3)), xmin = 0.104, xmax = 0.12, color='g')
    ax.hlines(y=func(0.115, *curve_fit(x, y3)), xmin = 0.104, xmax = 0.115, color='g')
    ax.hlines(y=func(0.12, *curve_fit(x, y2)), xmin = 0.104, xmax = 0.12, color = 'b')
    ax.hlines(y=func(0.115, *curve_fit(x, y2)), xmin = 0.104, xmax = 0.115, color='b')
    plt.text(0.123, 40.5, r'${\rm Z/ \gamma^\star + \geq 3\,jets}$', rotation =25)
    plt.text(0.123, 30, r'${\rm Z/ \gamma^\star + \geq 2\,jets}$', rotation=18)
    gap1 = round(func(0.12, *curve_fit(x, y3)) - func(0.115, *curve_fit(x, y3)), 2)
    gap2 = round(func(0.12, *curve_fit(x, y2)) - func(0.115, *curve_fit(x, y2)), 2)
    plt.text(0.105, (func(0.12, *curve_fit(x, y3))+func(0.115, *curve_fit(x, y3)))/2, str(gap1))
    plt.text(0.105,(func(0.12, *curve_fit(x, y2))+func(0.115, *curve_fit(x, y2)))/2, str(gap2))
    plt.annotate(s='', xy=(0.1045,func(0.12, *curve_fit(x, y2))), xytext=(0.1045,func(0.115, *curve_fit(x, y2))), arrowprops=dict(arrowstyle='<->', color='r'))
    plt.annotate(s='', xy=(0.1045,func(0.12, *curve_fit(x, y3))), xytext=(0.1045,func(0.115, *curve_fit(x, y3))), arrowprops=dict(arrowstyle='<->', color='r'))
    plt.annotate(s='', xy=(0.115,20.5), xytext=(0.12,20.5),arrowprops=dict(arrowstyle='<->', color='r') )
    plt.text(0.1165,21, '0.005')
    plt.xlabel(r'${\rm \alpha _s}$')
    plt.ylabel( r'$\rm Total \ \rm cross \ \rm section, \ \sigma,\ (pb)}$')

    
    #r'$\rm p_{\rm T}^{\rm jet} \,(leading\, jet) \,[GeV]}$
    plt.show()

def func(x,a, b):
    return a*x+b


def curve_fit(x,y):

    params, params_covariance= optimize.curve_fit(func,x,y)
    #print params
    return params 


def covst(e):
    return np.diag(e**2)

def cov_other(e,y):
    C = np.zeros((len(y),len(y)), dtype=np.float64)
    for i in range(len(y)):
        for j in range(len(y)):
            C[i,j] = e[i]*e[j]
    return C


def total_errorst(e, bin):
    return np.matmul(bin,np.matmul(covst(e),bin.T)) 

def total_error_other(e, bin,y):
    return np.matmul(bin,np.matmul(cov_other(e,y),bin.T)) 

'''The code below is used to calculate the C_exp matirix for the electron channel and stored
in Cov_all with the keys written below, want to compare the diagonals not squared with no luminosity included'''
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
            X[i,j] = (Cov_syst(e)[i,j])/100 * np.sqrt(y[i]*y[j])
    return X

if __name__ =="__main__": 
    y2 = np.array(data_z2['Cross_Section_Zll [pb]'])
    y3 = np.array(data_z3['Cross_Section_Zll [pb]'])
    y4 = np.array(data_z4['Cross_Section_Zll [pb]'])
    #statistical error
    st2 = np.array(data_z2['stat +'])
    st3 = np.array(data_z3['stat +'])
    st4 = np.array(data_z4['stat +'])
    #systematic error
    sy2 = np.array(data_z2['syst +'])
    sy3 = np.array(data_z3['syst +'])
    sy4 = np.array(data_z4['syst +'])
     #lumi error
    l2 = np.array(data_z2['lumi +'])
    l3 = np.array(data_z3['lumi +'])
    l4 = np.array(data_z4['lumi +'])
    binx= np.array([20,20,20,20,20,20,20,30,50,50,50,50,100])
    binx4= np.array([20,20,20,20,20,20,20,30,50,50,100,100])



    #full systematic error breakdown. 
    er_2e = np.array(error_2es.iloc[1:44,1:])
    lum = np.array(error_2es.iloc[44:45,1:])
    their_er2 = np.array(data_e2['syst +'])
    their_lumi = np.array(data_e2['lumi +'])
    their_Data = np.array(data_e2['Cross_Section_Zee [pb]'])
    
    ''''
    my_error = np.diag(Cov_syst100(er_2e,their_Data)) #this code finds the ratio of my electron errors to the ones quoted in hep.
    their_error =  their_er2
    print 'theirerror', their_error
    print 'myerror', my_error
    print  their_error/my_error
    my_lum = np.diag(Cov_syst100(lum,their_Data))
    print 'mylum', my_lum
    #print their_lumi/my_lum #now gives much better results. '''

    '''This code gives the total cross section errors for comparison with the ones quoted in HEP '''
    #print np.sqrt(total_errorst(st3,binx))
    #print np.sqrt(total_error_other(sy3,binx,y3))
    #print np.sqrt(total_error_other(l3,binx,y3))
    #print tot_xsec_exp(y2,binx)
    #print tot_xsec_exp(y3,binx)
    #print tot_xsec_exp(y4,binx4)
    
    

    #print np.sqrt(np.sum(st2**2))

    
    #start = 1
    #stop = 22
    #step = 1
    #pdf = 'MSTW08cl0'
    
    #n = 22
    #scalevariation(binx, start, stop, step, pdf,n)

    x = np.arange(107,128)/1000.
   

    tot_Vs_as_plt(x)

Cov_all.close()
