'''This module evaluates the error vectors for NN pdfs'''
import numpy as np
import shelve
import pandas as pd

jet2 = pd.read_csv('HEPData-ins1514251-v2-Table_66.csv')
jet3 = pd.read_csv('HEPData-ins1514251-v2-Table_67.csv')


#def run_scale_variation(old_pdf_name, new_pdf_name, old_pdf_name1):
#def run_scale_variation(pdf, start, stop, step):
def run_scale_variation(label,y2,y3):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
               if a == 1 and b == 1:
                theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet')
                theory_z3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet')
                PDF_er2 = shelve.open('PDF_er2('+str(a)+','+str(b)+')')
                PDF_er3 = shelve.open('PDF_er3('+str(a)+','+str(b)+')')
                #shelve_func(theory_z2, old_pdf_name1, y2, PDF_er2, new_pdf_name)
                #shelve_func(theory_z3, old_pdf_name1, y3, PDF_er3, new_pdf_name)

                PDF_er2.close()
                PDF_er3.close()

                print a,b
               else:
                print a,b
                theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_('+str(a)+',' +str(b)+ ')')
                theory_z3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet_('+str(a)+',' +str(b)+ ')')

                yNN = shelve.open('yNN_('+str(a)+',' +str(b)+ ')')
                #PDF_er2 = shelve.open('PDF_er2('+str(a)+','+str(b)+')')
                #PDF_er3 = shelve.open('PDF_er3('+str(a)+','+str(b)+')')
                for i in np.array([115,117,118,119,121]):
                    values  = Av(theory_z2, label, y2, i)
                    yNN[label+str(i)+'_j2'] = values
                for i in np.array([115,117,118,119,121]):
                    values = Av(theory_z3, label,y3, i)
                    yNN[label+str(i)+'_j3'] = values

                #shelve_func_pc(corr2, PDF_er2, pdf, start, stop, step)
                #shelve_func_pc(corr3, PDF_er3, pdf, start, stop, step)

                #shelve_func(theory_z2, old_pdf_name, y2, PDF_er2, new_pdf_name)
                #shelve_func(theory_z3, old_pdf_name, y3, PDF_er3, new_pdf_name)
                yNN.close()
                theory_z2.close()
                theory_z3.close()
                #theory_z4.close()

                #PDF_er2.close()
                #PDF_er3.close()
                #Cpd4.close()


def si(jet, label, y, k, i): #label = 'NNPDF23_0'
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    # k is the index of the alpha value, i is the index of the index of the pdf set
    return np.array(jet[label + str(k) + '_' + str(i)+''])[1:y+1] #for particular pdf parameter change k, is a vector length y

def Av(jet,label,y,k): #finds the average value for the statistical pdf's'
    Av = np.zeros((y), dtype=np.float64)
    for i in range(101):
        Av = si(jet, label, y, k, i) + Av
    #print Av
    return Av/100 #divide by n

def yNNc(jet,label,y, j):
    
    for i in np.array([115,117, 118, 119, 121]):
        values = Av(jet, label, y, i)
        yNN[label+str(i) + j] = values
    

def Diff(jet,label, y, k, i ):
    return (si(jet, label, y, k, i) - Av(jet, label, y, k))**2

def STDdev(jet, label, y, k): #the errors are the standard deviation, and combined them in the same way as the hessian. 
    s = np.zeros((y), dtype=np.float64)
    for i in range(101):
        s = Diff(jet, label, y, k, i) + s
    return np.sqrt(s/99) #divide by n-1, the -1 is negligible for large n though.

def Cpdf(jet,label,y,k):
    Cpdf = np.zeros((y,y), dtype = np.float64)
    for i in range(y):
        for j in range(y):
            Cpdf[i,j] = STDdev(jet,label,y,k)[i] * STDdev(jet,label,y,k)[j] #outer product of the errors to find the Cpdf
    return Cpdf

def shelve_func(jet, label, y,PDF_er,pdf): #store the Cpdf matrices in a shelve.

    for k in range(114,125): #(114-125 for NN23, 115, 117, 118, 119, 121 for NN30
        values = STDdev(jet,label,y,k)
        print np.shape(values)
        PDF_er[pdf+ str(k)+''] = values  #returns the errors for all keys and as values. 

def pc_errors(corr, PDF_er, pdf,k):
    return PDF_er[pdf+str(k)+''] * corr

def shelve_func_pc(corr, PDF_er, pdf, start, stop, step):
    for k in range(start, stop, step):
        values = pc_errors(corr, PDF_er, pdf, k)
        PDF_er[pdf+str(k)+''] = values

if __name__ =="__main__":
    y2 = 13 
    y3 = 13
    y4 = 12  

    corr2 = np.array(jet2['NP_corrections_Zll [none]'])
    corr3 = np.array(jet3['NP_corrections_Zll [none]'])

    old_pdf_name = 'NNPDF23_0'
    new_pdf_name = 'NNPDF23_0' 
    old_pdf_name1 = 'NNPDF23_0'

    label = 'NNPDF30_0' 
    start = 118
    stop = 119
    step = 1


    theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_(0.5,0.5)')
    k = 115           
    Av(theory_z2,label,y2,k)
    theory_z2.close()
    
    
    run_scale_variation(label,y2,y3)
    #run_scale_variation(pdf, start, stop, step)
    #run_scale_variation(old_pdf_name, new_pdf_name, old_pdf_name1)
    #yNNc(theory_z2, 'NNPDF30_0', y2, '_j2')
    #yNNc['NN30_0118_j2'] = Av(theory_z2, 'NNPDF30_0', y2, 118) #NEED TO CALCULATE 


