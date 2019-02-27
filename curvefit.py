'''this module peforms a curve fitting to a chi squared funtion to find the best va;ue of alpha s. 
there will be an error associated with this fit and hence the parameters used. f(a(ur,uf)), 
minimise this f to find as(ur,uf). Also plots all the scales on one plot.'''

import numpy as np
import scipy
from scipy import optimize
from scipy.optimize import minimize
import shelve
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import itertools






def test_func(x,a,b,c,d):
    return a*x**3 +b*x**2 +c*x + d


def curve_fit(x,y):

    params, params_covariance= optimize.curve_fit(test_func,x,y)
    #print params
    return params 


def sort(ax,y,yth, ykey, jetst, Cpdjet, pdf, Ctj):

    lists = sorted(itertools.izip(*[ax, CHI_SQ(y, yth, ykey, jetst, Cpdjet, pdf, Ctj)]))
    return np.array(list(itertools.izip(*lists)))

def CHI_min(pdf, a_s, n):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    X =  np.zeros((3,3), dtype= np.float64)
   
    for i in range(3):
        for j in range(3):
            a = [0.5,1,2]
            b = [0.5,1,2]

            
            CHI_square = shelve.open('CHI_Squared_('+str(a[i])+','+str(b[j])+')')
            params = curve_fit(a_s, CHI_square[pdf][0:n])
            a_smin = xmin(*params)
            X[i,j] = a_smin
    
    
            Chi_min = test_func(a_smin, *curve_fit(a_s, CHI_square[pdf][0:n]))
            plt.scatter(a_smin, Chi_min, marker = 's')
            plt.plot(a_s, test_func(a_s, *curve_fit(a_s, CHI_square[pdf][0:n])), label='('+str(a[i])+','+str(b[j])+')')
            plt.scatter(a_s, CHI_square[pdf][0:n], marker='x', color='k')
    
    
    legend = ax.legend(loc='upper left', fontsize='small')
    plt.xlabel(r'${\rm \alpha _s}$')
    plt.ylabel(r'${\chi ^2(\alpha_s)}$')
    print 'alpha_s for ' + pdf + ' is', X[1,1]
    print 'alpha_upper_error is ', X[2,2] - X[1,1]
    print 'alpha_lower_error is', X[0,0] - X[1,1]
    plt.xlim(0.110,0.130)
    plt.savefig(pdf+'_curvefit.pdf')
    plt.show()     
    CHI_square.close()

def scale_error_uncorrelated(pdf, n, a_s):
    A_s = shelve.open('Alpha_S')  
    a = [0.5,1,2]
    b = [0.5,1,2]
    X =  np.zeros((3,3), dtype= np.float64)
    for i in range(3):
        for j in range(3):
            if (i==0 and j==2) or (i==2 and j==0):
                a,b
            
            else:
                CHI_square = shelve.open('CHI_Squared_('+str(a[i])+','+str(b[j])+')')
                params = curve_fit(a_s, CHI_square[pdf][0:n])
                a_smin = xmin(*params)
                X[i,j] = a_smin
                Chi_min = test_func(a_smin, *curve_fit(a_s, CHI_square[pdf][0:n]))

    print X
    #print np.nanmax(X) - X[1,1]
    #print X[1,1] - np.nanmin(X[np.nonzero(X)])
    #A_s[pdf] = np.append(A_s[pdf], [np.nanmax(X) - X[1,1], X[1,1] - np.nanmin(X[np.nonzero(X)])])
    #print A_s[pdf]
    #print np.shape(A_s[pdf]) 
    

def scale_error_correlated(pdf, n, a_s):
    A_s = shelve.open('Alpha_S')  
    a = [0.5,1,2]
    b = [0.5,1,2]
    X =  np.zeros((3,3), dtype= np.float64)
    for i in range(3):
        for j in range(3):
            if (i==j):
                CHI_square = shelve.open('CHI_Squared_('+str(a[i])+','+str(b[j])+')')
                params = curve_fit(a_s, CHI_square[pdf][0:n])
                a_smin = xmin(*params)
                X[i,j] = a_smin
                Chi_min = test_func(a_smin, *curve_fit(a_s, CHI_square[pdf][0:n]))
            
            else:
                a,b
    print X
    print np.nanmax(X) - X[1,1]
    print X[1,1] - np.min(X[np.nonzero(X)])
    A_s[pdf] = np.append(A_s[pdf], [np.nanmax(X) - X[1,1], X[1,1] - np.nanmin(X[np.nonzero(X)])])
    print A_s[pdf]
    print np.shape(A_s[pdf])
    A_s.close() 

def Chi_error_and_Result(pdf, a_s, n):

    A_s = shelve.open('Alpha_S')
    a = 1
    b = 1
    X = np.zeros((3))
    CHI_square = shelve.open('CHI_Squared_('+str(a)+','+str(b)+')')
    params = curve_fit(a_s, CHI_square[pdf][0:n])
    a_smin = xmin(*params)
    Chi_min = test_func(a_smin, *curve_fit(a_s, CHI_square[pdf][0:n]))

    X[0] = a_smin
    Chi_1sigma = Chi_min + 1
    print Chi_min
    print Chi_1sigma
    print a_smin
    # one_sigma = np.roots(test_func(a_s, *curve_fit(a_s, CHI_square[pdf][0:n])) - Chi_1sigma)
    one_sigma = np.roots([params[0],params[1], params[2], params[3] - Chi_1sigma])
    X[1] = one_sigma[1] - a_smin
    X[2] = a_smin - one_sigma[2]
    print one_sigma
    print X #X is alphas upper 1 sigma error, lower 1 sigma error.
    #A_s[pdf] = X

    A_s.close()


def xmin(a,b,c,d):

    return (-b + np.sqrt(b**2 -3*a*c))/(3*a)

def Chi_source_min(pdf, a_s, n, type_s):
    a=1
    b=1
    A_s = shelve.open('Alpha_S') #loading the alpha_s DATA
    Chi_sour = shelve.open('Chi_source') #loading the chi squared tests for each error source
    
    CHI_square = shelve.open('CHI_Squared_('+str(a)+','+str(b)+')')
    params = curve_fit(a_s, Chi_sour[type_s+pdf]) #parameters for  chi squared test.
    a_smin = A_s[pdf][0]
    Chi_min_sour = test_func(a_smin, *curve_fit(a_s, Chi_sour[type_s+pdf]))
    Chi_min = test_func(a_smin, *curve_fit(a_s, CHI_square[pdf][0:n]))
    error = Chi_min_sour/Chi_min
    X = np.zeros((2), dtype=np.float64)
    print type(A_s[pdf])
    #A_s[pdf] = np.append(A_s[pdf], [np.sqrt(error * A_s[pdf][1]**2), np.sqrt(error * (A_s[pdf][2])**2)])
    
    print A_s[pdf]

    A_s.close()

    






if __name__ =="__main__": 

 pdf = 'NNPDF30_0'
 #a_s = np.arange(114,124)/1000.
 
 a_s =  np.array([0.115,0.117,0.118,0.119,0.121])
 n = 5
 #type_s = 'theory_' #need to do theory now
 #print np.shape(a_s)
 #Chi_error_and_Result(pdf, a_s, n)
 #CHI_min(pdf, a_s, n)
 #scale_error_uncorrelated(pdf, n, a_s)
 #scale_error_correlated(pdf, n, a_s)
 #Chi_source_min(pdf, a_s, n, type_s) #maybe do as a percentage of the total uncertainty
 #scale_variation(pdf, a_s)

 A_s = shelve.open('Alpha_S')

 #arr = A_s[pdf]
 #x = np.delete(arr, [11,12], None)
 #print x
 #A_s[pdf] = np.delete(arr, [11,12], None)
 #print A_s[pdf]

 #A_s.close()
 
 A_s['text'] = 'alpha_s', 'total_upper_error', 'total_lower_error', 'PDF_upper_error', 'PDF_lower_error', 'Exp_upper_error', 'Exp_lower_error', 'theory_ue', 'theory_le', 'scale_uncor_ue', 'scale_uncorr_le', 'scale_corr_ue', 'scale_corr_le'
 
 print  A_s['text']

