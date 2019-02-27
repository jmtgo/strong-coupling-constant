'''This module evaluates the error vectors of the pdf's from the individual error matrices'''

import shelve
import numpy as np

def run_scale_variation(pdf, start,stop, step):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:

                PDF_er2 = shelve.open('PDF_er2('+str(a)+','+str(b)+')')
                PDF_er3 = shelve.open('PDF_er3('+str(a)+','+str(b)+')')
                
                PDF_er = shelve.open('PDF_er('+str(a)+','+str(b)+')')
                
                for i in range(start, stop, step):
                    A = np.append(PDF_er2[pdf+str(i)], PDF_er3[pdf+str(i)])
                    B = np.append(A, PDF_er2[pdf+str(i)])
                    values = np.append(B, PDF_er3[pdf+str(i)])
               
                    PDF_er[pdf+str(i)] = values
                    print np.shape(values)




                PDF_er.close()
                PDF_er2.close()
                PDF_er3.close()
                #Cpd4.close()

   

if __name__ =="__main__": 

    y2 = 13
    y3 = 13
    
    pdf = 'CT10nlo0'

    start = 127
    stop = 128
    step = 1

    run_scale_variation(pdf, start,stop, step)

    
    
    
   