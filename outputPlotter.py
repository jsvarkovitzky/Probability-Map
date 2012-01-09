#This program reads in most2geoclaw.py output and plots it without any manipulation
#Jonathan Varkovitzky
#1-09-2012

from numpy import *
from matplotlib import *
from matplotlib import rc
rc("text",usetex=True)
from pylab import *
from scipy.special import *

###########################################
## Read in a single time step of one run ##
###########################################

def readIn(frameNum,run,n,m,k):
    print "Reading in frame %s of run %r" % str(frameNum,run)



##################
## Main Program ##
##################

readIn(1,2,10,10,10)
