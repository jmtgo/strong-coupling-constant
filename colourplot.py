'''This module evaluates a colour plot so as to show the substructre of the C_Exp matrix, the errors used
here are percentage error so as to display the information most clearly'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.colors as colors
import shelve


Cov_all = shelve.open('Cov_all')
yexp4 = Cov_all['FULL_xsec_vec']
C_Exp = Cov_all['FULL_ALLTHREE']



def re_scale(y):
    X = np.zeros((np.shape(C_Exp)), dtype=np.float64)
    for i in range(len(C_Exp[0,:])):
        for j in range(len(C_Exp[0,:])):
            X[i,j] = C_Exp[i,j]/(y[i]*y[j])
    return X  

def colourplot(x,y, array):

    plt.subplot(111)
    pcm = plt.pcolor(x,y,array, norm=colors.LogNorm(vmin = np.amin(array) ,vmax= np.amax(array)),cmap= 'Purples')
    cbar = plt.colorbar(pcm)
    cbar.set_label('value of element of Cov_matrix')
    plt.show()

if __name__ =="__main__": 
    
    x = np.arange(76)
    y=x
    print np.shape(C_Exp)
    array = re_scale(yexp4)

    colourplot(x,y,array)