'''This module evaluates the Cpdf covariance matrix for any pdf (not NN type), of the principle as value'''


import numpy as np
import shelve
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors


#a = #(1,0.5) (1,2) (0.5,0.5) (0.5,1) (0.5,2) (2,1) (2,2) (2,0.5)
#b =


#theory_z4 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/fourjet')
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
#data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')
#Cpd4 = shelve.open('Cpdf4')

def run_scale_variation(old_pdf_name, new_pdf_name, start, stop, pc_pred_as):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
            if a == 1 and b == 1:
                theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet')
                theory_z3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet')
                print a,b
                PDF_er2 = shelve.open('PDF_er2('+str(a)+','+str(b)+')')
                PDF_er3 = shelve.open('PDF_er3('+str(a)+','+str(b)+')')
                pc2 = shelve.open('pc2__('+str(a)+','+str(b)+')')
                pc3 = shelve.open('pc3__('+str(a)+','+str(b)+')')

                #print np.shape(Av_e(theory_z2, old_pdf_name1,y2))
                #PDF_er2[new_pdf_name] = Av_e(theory_z2, old_pdf_name1,y2)
                #PDF_er3[new_pdf_name] = Av_e(theory_z3, old_pdf_name1,y3)

                for i in range(start, stop):
                    PDF_er2[new_pdf_name + '0' + str(i)]= PDF_pc(theory_z2,old_pdf_name, PDF_er2, pc_pred_as, pc2, new_pdf_name, y2, i)
                    PDF_er3[new_pdf_name + '0' + str(i)]= PDF_pc(theory_z3,old_pdf_name, PDF_er3, pc_pred_as, pc3, new_pdf_name, y3, i)

                PDF_er2.close()
                PDF_er3.close()
                pc2.close()
                pc3.close()

            else:
                print a,b
                theory_z2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_('+str(a)+',' +str(b)+ ')')
                theory_z3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet_('+str(a)+',' +str(b)+ ')')

                Cpd2 = shelve.open('Cpdf2_('+str(a)+',' +str(b)+ ')') #will need to specify the PDF later in the code
                Cpd3 = shelve.open('Cpdf3_('+str(a)+',' +str(b)+ ')') #(1,1) (1,0.5) (1,2) (0.5,0.5) (0.5,1) (0.5,2) (2,1) (2,2) (2,0.5)
                PDF_er2 = shelve.open('PDF_er2('+str(a)+','+str(b)+')')
                PDF_er3 = shelve.open('PDF_er3('+str(a)+','+str(b)+')')
                pc2 = shelve.open('pc2__('+str(a)+','+str(b)+')')
                pc3 = shelve.open('pc3__('+str(a)+','+str(b)+')')

                for i in range(start, stop):
                    PDF_er2[new_pdf_name +'0' + str(i)]= PDF_pc(theory_z2,old_pdf_name, PDF_er2, pc_pred_as, pc2, new_pdf_name, y2, i)
                    PDF_er3[new_pdf_name + '0' + str(i)]= PDF_pc(theory_z3,old_pdf_name, PDF_er3, pc_pred_as, pc3, new_pdf_name, y3, i)
                #PDF_er2[new_pdf_name] = Av_e(theory_z2, old_pdf_name,y2)
                #PDF_er3[new_pdf_name] = Av_e(theory_z3, old_pdf_name,y3)
                

                theory_z2.close()
                theory_z3.close()
                #theory_z4.close()
                Cpd2.close()
                Cpd3.close()
                PDF_er2.close()
                PDF_er3.close()
                pc2.close()
                pc3.close()
                #Cpd4.close()


def si(jet, label,y,k):
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet[label + str(k) + ''])[1:y+1] #for particular pdf parameter change k, is a vector length y

def Up_Diff(jet,label, y, k): #returns a vector length y.
    D = np.amax([np.zeros(y), (si(jet, label, y, k+1)- si(jet,label, y, 0)), (si(jet, label, y, k)- si(jet,label, y, 0))], axis=0)
    return D #finds the maximum of 0, si0-sik+, si0+sik-, where the k's come in pairs for plus and minus the parameter change. 

def UpE_sq(jet, label, y): #returns a vector
    UE = np.zeros(y)
    for k in range(1,41,2):
        
        '''' #53 for CT10, 57 for CT14, 41 for MSTW, 51 for MMHT'''


        UE = (Up_Diff(jet, label, y, k))**2 + UE
    return UE #square all the D's for each k pair and adds them together.

def Lo_Diff(jet,label, y, k): #returns a vector length y
    D = np.amin([np.zeros(y), (si(jet, label, y, k+1)+ si(jet,label, y, 0)), (si(jet, label, y, k)+ si(jet,label, y, 0))], axis=0)
    return D #finds the minimum of 0, si0+sik+, si0 + sik-

def LoE_sq(jet, label, y): #returns a vector length y
    LE = np.zeros(y)
    for k in range(1,41,2): 
        
        
        '''#53 for CT10, 57 for CT14, 41 for MSTW, 51 for MMHT'''


        LE = (Lo_Diff(jet, label, y, k))**2 + LE
    return LE #square all the D's for each k pair and adds them together.

def Av_e(jet, label, y): #averages the upper and lower errors. 
    return (np.sqrt(LoE_sq(jet,label,y)) + np.sqrt(UpE_sq(jet, label, y)))/2
      

def Cpdf(jet,label,y):
    Cpdf = np.zeros((y,y), dtype = np.float64)
    for i in range(y):
        for j in range(y):
            Cpdf[i,j] = Av_e(jet,label,y)[i] * Av_e(jet,label,y)[j] #outer product of the errors to find the Cpdf
    return Cpdf


def PDF_percent(shelve, key, pred_ao_shelve, pdf,y, i):
    return pred_ao_shelve[pdf]/shelve[key + str(i)][1:y+1]

def PDF_pc(shelve,key, pred_ao_shelve, pc_pred_as, pc_shelve, pdf,y, i):
    #shelve = theory_zi predictions from fastnlo directly
    #key = key for theory_zi
    #pred_ao_shelve = PDF_eri shelf with PDF errors in for alpha 0.
    #pdf = key for PDF_er shelf number i.
    # y  is length of array for slicing the fastnlo shelf.
    # pc_shelve = perturbative correction shelf
    # pc_pred_as = key for perturbative correction shelf. 
    return PDF_percent(shelve, key, pred_ao_shelve, pdf,y, i) * pc_shelve[pc_pred_as+str(i)]



if __name__ =="__main__":
    x2 = np.array(data_z2['pT(jet) [GeV]'])
    x3 = np.array(data_z3['pT(jet) [GeV]'])
    #x4 = np.array(data_z4['pT(jet) [GeV]'])
    y2 = 13
    y3 = 13
    y4 = 12   
   
    #old_pdf_name1 = 'MSTW0868cl'
    old_pdf_name = 'MMHT14_000' # key for fastnlo shelf
    new_pdf_name = 'MMHT14cl' #name of the PDF_0 shelf key
    start = 1
    stop = 21
    pc_pred_as = 'MMHT14cl0' #key for perturbative correction shelf.
    
    #PDF_pc(theory_z3,old_pdf_name, PDF_er3, pc_pred_as, pc3, new_pdf_name, y3, i)
    run_scale_variation(old_pdf_name, new_pdf_name, start, stop, pc_pred_as)

    #Cpd4['MSTW08cl'] = Cpdf(theory_z4,'MSTW0868cl', y4)
  
   

