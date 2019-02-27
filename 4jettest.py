'''This module compares the 4 jet cross section with the real data, to see if the right data is
being generated and the large error comes from a small data set'''

import numpy as np
import shelve
import matplotlib.pyplot as plt



def tot_xsec_exp(y,bin):
    return np.sum(y*bin)

xsec_pdf = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/fourjet')
x2sec_pdf = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/twojet')
x3sec_pdf = shelve.open('/Users/Georgia/Documents/DurhamAcademic/physics/linux-shelves/threejet')

xm = np.array([0.00294,0.00836,0.01013,0.01022,0.00831,0.0066,0.00547,0.00412,0.00274,0.00165,0.000724,0.000267]) #4 jets 
x2m = np.array([0.31,0.3351,0.2236,0.1473,0.0965,0.0639,0.0444,0.0296,0.01553,0.00735,0.0037,0.00217,0.000938])
x3m = np.array([0.0281,0.0557,0.0524,0.043,0.0315,0.0231,0.0172, 0.01273,0.00749,0.00407,0.00204,0.00122,0.000561])
binx= np.array([20,20,20,20,20,20,20,30,50,50,50,50,100])
binx4= np.array([20,20,20,20,20,20,20,30,50,50,100,100])

x = np.arange(0,len(xm))


xCT10 = np.array(xsec_pdf['0118'])[1:len(xm)+1] #4jet
x2CT10 = np.array(x2sec_pdf['0118'])[1:len(x2m)+1]
x3CT10 = np.array(x3sec_pdf['0118'])[1:len(x3m)+1]


#print len(xCT10)
fig = plt.figure()
ax = fig.add_subplot(111)
#plt.scatter(x,xm/xCT10)
#plt.scatter(x, xt/x10pdf, color='r')
#plt.scatter(x,xt/x30pdf, color='k')
#plt.show()

print tot_xsec_exp(xm,binx4)/tot_xsec_exp(xCT10, binx4) #gives the ratio of 4 jet measured total cross section to the total cross section of 4 jet from fastnlo predictions.
print tot_xsec_exp(x2m,binx)/tot_xsec_exp(x2CT10, binx)
print tot_xsec_exp(x3m,binx)/tot_xsec_exp(x3CT10, binx) #seems to be a factor of 8 out. 

xsec_pdf.close()
x2sec_pdf.close()
x3sec_pdf.close()