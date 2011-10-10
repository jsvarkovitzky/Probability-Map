#This code reads in the maximum wave height from provided files and combines them based on the probabilities of events to generate a probabilistic hazard map

#Jonathan Varkovitzky
#10-3-2011

from numpy import *
from matplotlib import *
from pylab import *
event1 = loadtxt('fort.fg01_0020',skiprows=9)


mu = 100
n = 200
m = 200

#Need to take the column of inundaiton data and convert it into a rectangular matrix
h = zeros((n,m))
for i in range (0,n-1):
    for j in range (0,m-1 ):
       h[i,j] = event1[i*n+j,3] #CHECK THESE INDICIES!!!
  

muMap = zeros((n,m))

for i in range(1,n):          #Loop through x direction
    for j in range(1,m):      #Loop through y direction
        if h[i,j] > 0.0:
            muMap[i,j] = muMap[i,j] + mu


matplotlib.pyplot.contour(h,100)
colorbar()
show()
