import numpy as np
import shelve
import pandas as pd
import matplotlib.pyplot as plt

data_z1 = pd.read_csv('HEPData-ins1514251-v2-Table_15.csv') #loading the data
data_z2 = pd.read_csv('HEPData-ins1514251-v2-Table_18.csv')
data_z3 = pd.read_csv('HEPData-ins1514251-v2-Table_21.csv')
data_z4 = pd.read_csv('HEPData-ins1514251-v2-Table_24.csv')

Covst = shelve.open('Cov_all')
Cpd2 = shelve.open('Cpdf2')
Cpd3 = shelve.open('Cpdf3')
Cpd4 = shelve.open('Cpdf4')
yNN = shelve.open('yNN_('+str(a)+',' +str(b)+ ')')


def run_scale_variation(old_pdf_name, new_pdf_name):
    '''Loops through the factorisation and renormalisation scales to create covariance matrices at each scale and store in a shelf'''
    for a in [0.5,1,2]:
        for b in [0.5,1,2]:
            if a == 1 and b == 1:
                # Do nothing in this case
                print a,b
            else:
                print a,b
                nocorr2 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet_('+str(a)+',' +str(b)+ ')')
                nocorr3 = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet_('+str(a)+',' +str(b)+ ')')

                Cpd2 = shelve.open('Cpdf2_('+str(a)+',' +str(b)+ ')') #will need to specify the PDF later in the code
                Cpd3 = shelve.open('Cpdf3_('+str(a)+',' +str(b)+ ')') #(1,1) (1,0.5) (1,2) (0.5,0.5) (0.5,1) (0.5,2) (2,1) (2,2) (2,0.5)
                

                theory_z2.close()
                theory_z3.close()
                #theory_z4.close()
                Cpd2.close()
                Cpd3.close()
                #Cpd4.close()
def Cov(i, jetst, Cpdjet, pdf):
    a_s = np.array([115, 117, 118, 119, 121])
    return np.array(Covst[jetst]) + np.array(Cpdjet[pdf + str(a_s[i])+'']) #pdf key fpr the Cpdf shelve

def CHI_SQ2(y, ykey, jetst, Cpdjet, pdf):
    X = np.zeros((5), dtype = np.float64)    #the size here should reflect the number of as values for each pdf set.
    a_s = np.array([115, 117, 118, 119, 121])
    for i in range(5): #ykey is the key of twojet, three jet etc of the desired pdf. #yth is the jet number type theory_z2
        X[i] = np.inner((y - np.array(yNN[ykey+str(a_s[i])+ '_j2'])), (np.inner(np.linalg.inv(Cov(i, jetst, Cpdjet, pdf)),(y - np.array(yNN[ykey+str(a_s[i])+'_j2'])).T)))
       
    return X

def CHI_SQ3(y, ykey, jetst, Cpdjet, pdf):
    X = np.zeros((5), dtype = np.float64)    #the size here should reflect the number of as values for each pdf set.
    a_s = np.array([115, 117, 118, 119, 121])
    for i in range(5): #ykey is the key of twojet, three jet etc of the desired pdf. #yth is the jet number type theory_z2
        X[i] = np.inner((y - np.array(yNN[ykey+str(a_s[i])+ '_j3'])), (np.inner(np.linalg.inv(Cov(i, jetst, Cpdjet, pdf)),(y - np.array(yNN[ykey+str(a_s[i])+'_j3'])).T)))
       
    return X


def CHI_graph(ax, ykey, pdf, y2, jet2st, Cpdjet2, y3, jet3st, Cpdjet3):
    fig1 =plt.figure(1)

    plt.scatter(ax,CHI_SQ2(y2, ykey, jet2st, Cpdjet2, pdf), marker = 'o', color = 'k')
    plt.scatter(ax,CHI_SQ3(y3, ykey, jet3st, Cpdjet3, pdf), marker = 'o', color = 'k')
    #plt.scatter(ax,(CHI_SQ('fourj', y4, yth4, C4)), marker = 'o', color = 'k')
    plt.plot(ax,CHI_SQ2(y2, ykey, jet2st, Cpdjet2, pdf), label = r'${\rm 2 \ jets}$' , color = 'b')
    plt.plot(ax,CHI_SQ3(y3, ykey, jet3st, Cpdjet3, pdf), label = r'${\rm 3 \ jets}$', color='r')
    #plt.plot(ax,CHI_SQ('fourj', y4, yth4, C4), label = r'${\rm 4 \ jets}$', color='g')
    plt.xlabel(r'${\rm \alpha _s}$')
    plt.ylabel(r'${\chi ^2(\alpha_s)}$')
    plt.legend(loc = 'upper left')
    plt.savefig(pdf + 'X^2')
    plt.show()


def si(jet_number_pdf, pdf,y,i,k):
   
    #slices the data as there are more theoretical values than measured values and want to take the correct corresponding bins
    return np.array(jet_number_pdf[pdf+str(i) + '_'+str(k)+''])[1:y+1]

def yNN_average2(pdf):
    values = np.zeros(13)
    for i in range(121,122):
        for j in range(13):
            for k in range(101):
                values = np.average(si(nocorr2, pdf,13,i,k)[j])
                yNN['NNPDF30_0' +str(i)+'_j2'] = values

def yNN_average3(pdf):
    values = np.zeros(13)
    for i in range(121,122):
        for j in range(13):
            for k in range(101):
                values[j] = np.average(si(nocorr3, pdf,13,i,k)[j])
                yNN['NNPDF30_0' +str(i)+'_j3'] = values


if __name__ =="__main__": 
    y1 = np.array(data_z1['Cross_Section_Zll [pb]'])
    y2 = np.array(data_z2['Cross_Section_Zll [pb]'])
    y3 = np.array(data_z3['Cross_Section_Zll [pb]'])
    y4 = np.array(data_z4['Cross_Section_Zll [pb]'])
    ax = np.array([115, 117, 118, 119, 121])


    #print np.shape(nocorr2['NNPDF23_0117_100'])
    yNN_average2('NNPDF30_0')
    yNN_average3('NNPDF30_0')
    #CHI_graph(ax, 'NNPDF30_0', 'NN30_0', y2, 'twoj', Cpd2, y3, 'threej', Cpd3)
    #print np.shape(yNN['NNPDF23_0114_j2'])
    #print np.shape(y2)

Covst.close()
Cpd2.close()
Cpd3.close()
Cpd4.close()
yNN.close()