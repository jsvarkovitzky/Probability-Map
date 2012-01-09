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

def readIn(frameNum,run,nx,ny):
    print "Reading in frame %s of run %r..." %(frameNum,run)
    outdir = 'run%s/_output/' %str(run).zfill(2)
    fname = outdir + 'fort.q%s' % str(frameNum).zfill(4)
    varname = 'fort_%s' % str(frameNum).zfill(4)
    exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))
    print "File imported sucessfully..."
    
    #Extract h and reshape into a matrix
    print "Extracting and reshaping h..."
    exec('y = %s[:,0].copy()'%varname)
    h = y.reshape(nx,ny,order='F')
    #Extract eta and reshape into a matrix
    print "Extracting and reshaping eta..."
    exec('z = %s[:,3].copy()'%varname)
    eta = z.reshape(nx,ny,order='F')
     
    print "Returning h and eta..."
    return(h,eta)

#########################
## Plotting Algorithms ##
#########################
def plotting(h,eta,frameNum,run):
    
    figure(1)
    clf()
    pcolormesh(h)
    colorbar()
    title("H of frame %s, run %r" %(frameNum,run))
    savefig("H_frame%s_run%r" %(frameNum,run))

    figure(2)
    clf()
    pcolormesh(eta)
    colorbar()
    title("Eta of frame %s, run %r" %(frameNum,run))
    savefig("Eta_frame%s_run%r" %(frameNum,run))
##################
## Main Program ##
##################

frameNum = 1
runNum = 17
nx = 271
ny = 192
k = 1

(h,eta) = readIn(frameNum,runNum,nx,ny)
plotting(h,eta,frameNum,runNum)
