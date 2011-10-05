#This program generates pseudo-maximum wave heights from a simulation to be used by probMaker.py

#Jonathan Varkovitzky
#10-2-2011

from scipy import *
from pylab import *
from matplotlib import *

plot = 0
n = 100
m = 100
k = 5
A = ((1,2,3,4))
#sigma = n/4
xCenter = n/2
yCenter = n/2

gauss = zeros((n,m,k))
for k in range(1,5):
    
    for i in range (0,n):
        for j in range (0,m):
            #        gauss_k[i,j] = A*exp(-((i-xCenter)**2+(j-yCenter)**2)/sigma)
            gauss[i,j,k] = A[k-1]*sin(i*0.01)*sin(j*.01)
            if plot == 1:
                matplotlib.pyplot.contour(gauss,100)
                show()
                

                

