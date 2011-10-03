from scipy import *
from pylab import *
from matplotlib import *

plot = 1
n = 100
m = 100
A = 10
#sigma = n/4
xCenter = n/2
yCenter = n/2
gauss1 = zeros((n,m))
for i in range (0,n):
    for j in range (0,m):
#        gauss1[i,j] = A*exp(-((i-xCenter)**2+(j-yCenter)**2)/sigma)
        gauss1[i,j] = A*sin(i*0.01)*sin(j*.01)
if plot == 1:
    matplotlib.pyplot.contour(gauss1,100)
    show()

#Write array to file
savetxt('runup1.txt',gauss1,fmt="%12.6G")
