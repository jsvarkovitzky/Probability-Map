#This code reads in the maximum wave height from provided files and combines them based on the probabilities of events to generate a probabilistic hazard map

#Jonathan Varkovitzky
#10-3-2011

from numpy import *
from matplotlib import *
from pylab import *
event1 = loadtxt('runup1.txt')
event2 = loadtxt('runup2.txt')
event3 = loadtxt('runup3.txt')
event4 = loadtxt('runup4.txt')

mu = (76,32,88,120)
n = 100
m = 100

eventk = zeros((n,m,4))
eventk[:,:,1] = event1[:,:]
eventk[:,:,2] = event2[:,:]
eventk[:,:,3] = event3[:,:]
eventk[:,:,4] = event4[:,:]

muMap = zeros((n,m))

for k in range(0,4):
    for i in range(0,n):
        for j in range(0,m):
            if eventk[i,j,k] > 1:
                muMap[i,j] = muMap[i,j] + mu[k]


matplotlib.pyplot.contour(muMap,100)
show()
