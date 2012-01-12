#This program reads in most2geoclaw.py output and plots it without any manipulation
#Jonathan Varkovitzky
#1-09-2012

from numpy import *
from matplotlib import *
from matplotlib import rc
rc("text",usetex=True)
from pylab import *
from scipy.special import *
import time
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

    print "Setting no data values to 10"
    h = where(h>-10**10,h,-10)
    eta = where(eta>-10**10,eta,-10)
     
    print "Returning h and eta..."
    return(h,eta)

####################
## mu Calculation ##
####################
def muCalc(zeta,zeta_i,nx,ny,nu,T):
    print "Computing mu"
    
    #Constant values obtained from Mofjeld et al. 2007
    MSL = 1.130
    MHHW = 2.095
    alpha = 0.170
    beta = 0.858
    C = 1.044
    CP = 0.707 
    alphaP = 0.056
    betaP = 1.119
    sigma_0 = 0.638 
    mu = zeros((nx,ny))
    

    for ix in range(0,nx):     #loop through x
        status(ix,nx)
        for iy in range(0,ny): #loop through y
                #The following equations are from Gonzalez et al 2009
            zeta_0j = zeta[ix,iy] + MSL + C*(MHHW-MSL)*exp(-alpha*(zeta[ix,iy]/sigma_0)**beta) #eqn 2b
            sigma_j = sigma_0-CP*sigma_0*exp(-1*alphaP*(zeta[ix,iy]/sigma_0)**betaP)          #eqn 2c
            mu[ix,iy] = mu[ix,iy] + 1./2*nu*(1-erf((zeta_i-zeta_0j)/(sqrt(2)*sigma_j)))   #eqn 3b/5
  
    return(mu)

###################
## P Calculation ##
###################

def pCalc(mu,nx,ny,T):
    print "Computing P(zeta>zeta_i)"
    P = zeros((nx,ny))
    P = ones((nx,ny))-exp(-1*mu*T) #eqn 6
    return(P)

#########################
## Plotting Algorithms ##
#########################
def plotting(h,eta,max_h,max_eta,frameNum,run):

    outdir = "_plots/"
    figure(1)
    clf()
    pcolormesh(h)
    colorbar()
    title("H of frame %s, run %r" %(frameNum,run))
    savefig(outdir + "H_frame%s_run%r" %(frameNum,run))

    figure(2)
    clf()
    pcolormesh(eta)
    colorbar()
    title("Eta of frame %s, run %r" %(frameNum,run))
    savefig(outdir + "Eta_frame%s_run%r" %(frameNum,run))

    figure(3)
    clf()
    pcolormesh(max_h)
    colorbar()
    title("Max H at frame %s, run %r" %(frameNum,run))
    savefig(outdir + "Max_H_frame%s_run%r" %(frameNum,run))

    figure(4)
    clf()
    pcolormesh(max_eta)
    colorbar()
    title("Max Eta at frame %s, run %r" %(frameNum,run))
    savefig(outdir + "Max_Eta_frame%s_run%r" %(frameNum,run))

#########################
## Plotting Algorithms ##
#########################
def probPlotter(max_h,mu,P):
    figure(5)
    clf()
    pcolormesh(mu)
    colorbar()
    title("mu")
    
    figure(6)
    clf()
    pcolormesh(P)
    colorbar()
    title("P")
    
    figure(7)
    clf()
    pcolormesh(max_h)
    colorbar()
    title("max_h")
    show()
################
## Status Bar ##
################

def status(n,N):
    n = float(n)
    N = float(N)
    percent = n/N*100
    sys.stdout.write("[==> ]%3d%%\r" %percent)
    sys.stdout.flush()


##################
## Main Program ##
##################
close('all')

frameNum = 59
runNum = 15
nx = 271
ny = 192
k = 1
T = 1
T_M = 520
nu = 520
zeta_i = 2.1
max_h = zeros((nx,ny))
max_eta = zeros((nx,ny))
mu = zeros((nx,ny))

for frameNum in range(1,60):
    (h,eta) = readIn(frameNum,runNum,nx,ny)
    max_h = where(max_h<h,h,max_h)
    max_eta = where(max_eta<eta,eta,max_eta)
    if frameNum == 35:
        plotting(h,eta,max_h,max_eta,frameNum,runNum)
 
muHNum = muCalc(max_h,zeta_i,nx,ny,nu,T)   
P = pCalc(muHNum,nx,ny,T)
probPlotter(max_h,muHNum,P)
