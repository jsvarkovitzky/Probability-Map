# This program uses functions defined in inundationMapper.py to generate
# a probability map of inundation.

# Jonathan Varkovitzky
# 12-31-2011

from numpy import *
from matplotlib import *
from matplotlib import rc
rc("text",usetex=True)
from pylab import *
from scipy.special import *
from inundationMapper import *

close ('all')
#test = 0 implies to run code on full set of  inundations
#test = 1 implies to run a test on a 2x2 identity matrix instead of inundation

test = 1
if test == 1:
    zeta_i = 2.1                       #meters of inundation
    fieldType = 1                      #1 = far-field, 2 = near-field
    ((nx,ny,runs)) = ((2,2,1))
    field_1 = array([[1,0],[0,1]])*100
    field_2 = array([[1,0],[0,1]])*100

    T_M = 520
    T = 1
    nu = ones((1,runs))*1./T_M

(P,mu_1) = tidalUncert(field_1,zeta_i,fieldType,nx,ny,runs,nu,T)
(P,mu_2) = tidalUncert(field_2,zeta_i,fieldType,nx,ny,runs,nu,T)

mu = mu_1 + mu_2
P = ones((nx,ny))-exp(-1*mu*T)
