'''This module creates the MEGA MATRIX used in the FINAL EXTRACTION, initially it will only 
consider 2 and 3 jets for the electron channel but will eventually incorporate all three numbers
of jets and both the electron and muon channel'''

import numpy as np
import shelve
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.colors as colors


#y data vector for all jets and channels.
Cov_all = shelve.open('Cov_all')
yexp4 = Cov_all['FULL_xsec_vec']

yexp = Cov_all['FULL_xsec_vec_no4']

#Experimental Cov matrix for each jet - no scale dependance

C_Exp = Cov_all['FULL_ALLTHREE_no4']



#Theoretical Cov matrix for each jet - no scale dependance
Ct = shelve.open('Ct')
theory2 = Ct['Ct2']
theory3 = Ct['Ct3']


#PDF predicitions and Covariance matirx
def run_scale_variation(pdf, theory2, theory3, M2, M3, M4, M5, yexp, start, stop, step, n):
#def run_scale_variation(start,stop,step,pdf, pc_pdf,y):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
            #yt2 = shelve.open('pc2__('+str(a)+','+str(b)+')') #theoretical predictions
            #yt3 = shelve.open('pc3__('+str(a)+','+str(b)+')')

            ypred = shelve.open('ypred_('+str(a)+','+str(b)+')')
            

            #for i in range(start,stop,step):
                #values = full_pred(i, yt2, yt3, pc_pdf)
                #print values
                #ypred[pdf+str(i)] = values
            
            
            #yt2.close()
            #yt3.close()

            #PDF_er = shelve.open('PDF_er('+str(a)+','+str(b)+')')
            Cpdf_large = shelve.open('Cpdf_large('+str(a) + ','+str(b)+')')
            #for i in range(start,stop,step):
            #     values2 = full_CPDF(pdf, PDF_er,i,y)
            #    Cpdf_large[pdf+str(i)] = values2

            #PDF_er.close()
            CHI_square = shelve.open('CHI_Squared_('+str(a)+','+str(b)+')')

            values = CHI_squared(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred, start, stop, step, n)
            CHI_square[pdf] = values
            


            ypred.close()
            Cpdf_large.close()
            CHI_square.close()

def full_pred(i, yt2, yt3, pc_pdf):
    A = np.append([yt2[pc_pdf+str(i)]],[yt3[pc_pdf+str(i)]])
    B =  np.append(A, [yt2[pc_pdf+str(i)]] )
    C = np.append(B, [yt3[pc_pdf+str(i)]] )
    return C

def full_Cexp(name):
    return Cov_all[name]

def full_Ctheory(theory2,theory3, M2, M3, M4, M5): 
    M2big = np.kron(M2, theory2)   
    M3big = np.kron(M3,theory3)
    M4big = np.kron(M4, theory2)
    M5big = np.kron(M5, theory3)
    M = (M2big + M3big + M4big + M5big)
    return M



def full_CPDF(pdf, PDF_er, i, y):
    X = np.zeros((y,y))
    for k in range(y):
        for j in range(y):
            X[j,k] = PDF_er[pdf+str(i)][j]*PDF_er[pdf+str(i)][k]
    return X


def full_cov(pdf, i, Cpdf_large, theory2, theory3, M2, M3, M4, M5):
    return Cpdf_large[pdf+str(i)] + full_Ctheory(theory2,theory3, M2, M3, M4, M5) + C_Exp

def CHI_squared(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred, start, stop, step, n):
    X = np.zeros((n), dtype = np.float64)
    #for i in range(5):
    #    ref_val = np.array([115,117,118,119,121])
    #    X[i] = np.inner((yexp-ypred[pdf+str(ref_val[i])]), (np.inner(np.linalg.inv(full_cov(pdf, ref_val[i], Cpdf_large, theory2, theory3, M2, M3, M4, M5)), (yexp-ypred[pdf+str(ref_val[i])]).T )))
    for i in range(start,stop,step):     
        X[i-start] =  np.inner((yexp-ypred[pdf+str(i)]), (np.inner(np.linalg.inv(full_cov(pdf, i, Cpdf_large, theory2, theory3, M2, M3, M4, M5)), (yexp-ypred[pdf+str(i)]).T )))
    return X

def Chi_sourceNN30(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred):
    X = np.zeros((5), dtype = np.float64)
 
    for i in range(5):
        ref_val = np.array([115,117,118,119,121])
        diffy = (yexp-ypred[pdf+str(ref_val[i])]) 
        print np.shape(diffy)
        Cf = full_cov(pdf, ref_val[i], Cpdf_large, theory2, theory3, M2, M3, M4, M5)
        Cs = C_Exp
        X[i] =  np.inner(diffy, (np.inner(np.linalg.inv(Cf), (np.inner(Cs, (np.inner(np.linalg.inv(Cf), (diffy.T ))))))))
    return X

def Chi_source(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred, start, stop, step, n):
    X = np.zeros((n), dtype = np.float64)
    
    for i in range(start,stop,step):
        
        diffy = (yexp-ypred[pdf+str(i)]) 
        print np.shape(diffy)
        Cf = full_cov(pdf, i, Cpdf_large, theory2, theory3, M2, M3, M4, M5)
        Cs = Cpdf_large[pdf+str(i)] #C_Exp
        X[i-start] =  np.inner(diffy, (np.inner(np.linalg.inv(Cf), (np.inner(Cs, (np.inner(np.linalg.inv(Cf), (diffy.T ))))))))
    return X

def Chi_source_Save(pdf, theory2, theory3, M2, M3, M4, M5, yexp, start, stop, step, n):
    a = 1
    b = 1

    ypred = shelve.open('ypred_('+str(a)+','+str(b)+')')
    Cpdf_large = shelve.open('Cpdf_large('+str(a) + ','+str(b)+')')

    Chi_sour = shelve.open('Chi_source')
    values = Chi_source(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred, start, stop, step, n)
    #values = Chi_sourceNN30(pdf, Cpdf_large, theory2, theory3, M2, M3, M4, M5, yexp, ypred)
    Chi_sour['PDF_'+ pdf] = values

    ypred.close()
    Cpdf_large.close()

            

def CHI_plot(pdfname, x):
     fig = plt.figure()
     ax = fig.add_subplot(1,1,1)
     for a in [0.5,1,2]:
        for b in [0.5,1,2]:
     
            CHI_square = shelve.open('CHI_Squared_('+str(a)+','+str(b)+')')
            print np.shape(CHI_square[pdfname])
            plt.scatter(x, CHI_square[pdfname][0:21])
            plt.plot(x, CHI_square[pdfname][0:21], label='('+str(a)+','+str(b)+')')
            plt.xlim(0.105,0.130)
     legend = ax.legend(loc='upper left', fontsize='small')
     plt.text(0.110,105, 'MSTW08')
     plt.text(0.110,98, r'${\rm \ (a,b) = (\mu_F \mu_R)}$')
     

     plt.xlabel(r'${\rm \alpha _s}$')
     plt.ylabel(r'${\chi ^2(\alpha_s)}$')

     plt.savefig(pdfname+'.pdf')

     CHI_square.close()

            
def colourplot(x,y, array):

    plt.subplot(111)
    pcm = plt.pcolor(x,y,array, norm=colors.LogNorm(vmin = np.amin(array) ,vmax= np.amax(array)),cmap= 'Purples')
    cbar = plt.colorbar(pcm)
    plt.xlim(0,51)
    plt.ylim(0,51)
    plt.xticks(color='None')
    plt.yticks(color='None')
    cbar.set_label('Covariance matrix')
    plt.savefig('colourplot_C_full.png')
    plt.show()
    
    


if __name__ =="__main__": 
    

    #start = 114
    #stop = 124
    #step = 1
    #pdf = 'NNPDF23_0'
    #pc_pdf = 'NNPDF23_0'
    #y = 52

    M2 = np.diag([1,0,0,0])
    M3 = np.diag([0,1,0,0])
    M4 = np.diag([0,0,1,0])
    M5 = np.diag([0,0,0,1])
    
    #Cpdf_large = shelve.open('Cpdf_large(1,1)')
    #n = (stop - start)/step 
    #x = np.arange(52)
    #y = np.arange(52)
    #e = full_cov(pdf, 118, Cpdf_large, theory2, theory3, M2, M3, M4, M5)
    #array = np.zeros((52,52), dtype=np.float64)
    #for i in range(52):
    #    for j in range(52):

                #array[i,j] = e[i,j]/np.sqrt(e[i,i]*e[j,j])
   
    

    #array = C_Exp 
    #print np.shape(array)
    #colourplot(x,y, e) 

    #Cpdf_large.close()
    #run_scale_variation(pdf, theory2, theory3, M2, M3, M4, M5, yexp, start, stop, step, n)
    #full_Ctheory(theory2, theory3, M2, M3)

    #x = np.arange(107,128)/1000.
    #x = np.array([0.115,0.117,0.118,0.119,0.121])
    #print np.shape(x)
    start = 1
    stop = 21
    step = 1
    n = stop-start
    pdf = 'MMHT14cl0'
    #CHI_square = shelve.open('CHI_Squared_(1,1)')
    Chi_source_Save(pdf, theory2, theory3, M2, M3, M4, M5, yexp, start, stop, step, n)
    #print np.shape(CHI_square[pdfname])
    #print x
    
    #run_scale_variation(start,stop,step,pdf, pc_pdf,y)
    #CHI_plot(pdfname, x)
    #full_yexp(y2edata, y3edata)
    Cov_all.close()
    Ct.close()
