#This code reads in the maximum wave height from provided files and combines them based on the probabilities of events to generate a probabilistic hazard map

#Jonathan Varkovitzky
#10-3-2011

from numpy import *
from matplotlib import *
from pylab import *
event1 = loadtxt('fort.fg01_0020',skiprows=9)


mu = (76,32,88,120)
n = 200
m = 200

#Need to take the column of inundaiton data and convert it into a rectangular matrix
h = zeros((n,m))
for i in range (1,n):
    for j in range (1,m):
       h[i,j] = event1[i+j-1,4]
  

muMap = zeros((n,m))

for k in range(0,1):              #Loop through invidual events
    for i in range(0,n):          #Loop through x direction
        for j in range(0,m):      #Loop through y direction
            if eventk[i,j,k] > 1:
                muMap[i,j] = muMap[i,j] + mu[k]


matplotlib.pyplot.contour(muMap,100)
show()
