#This code reads in the maximum wave height from provided files and combines them based on the probabilities of events to generate a probabilistic hazard map

#Jonathan Varkovitzky
#10-3-2011

from numpy import *
from matplotlib import *
from pylab import *
event1 = loadtxt('fort.fg01_0007',skiprows=9)


mu = 100
n = 200
m = 200

#Need to take the column of inundaiton data and convert it into a rectangular matrix
h = zeros((n,m))
topo = zeros((n,m))
nanZone = zeros((n,m))
for i in range (0,n-1):
    for j in range (0,m-1 ):
        h[i,j] = event1[i*n+j,4]      #max wave height
        topo[i,j] = event1[i*n+j,3]   #Topography
        if event1[i*n+j,4] not in range(-10**3,10**3):
            nanZone[i,j] = 1
muMap = zeros((n,m))

for i in range(1,n):          #Loop through x direction
    for j in range(1,m):      #Loop through y direction
        if h[i,j] > 0.0:
            muMap[i,j] = muMap[i,j] + mu

figure(1)
matplotlib.pyplot.contour(h,100)
colorbar()


figure(2)
matplotlib.pyplot.contour(topo,100)
colorbar()

figure(3)
matplotlib.pyplot.contour(nanZone,100)
colorbar()


show()

