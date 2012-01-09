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

########################################################
## Iterate over each time step to find max inundation ##
########################################################

def maxH( nfiles,run,n,m,k):
    print 'calculating maxH for run '+str(run)
    fnums = range(1,nfiles)
    max_h = zeros((n,m))
    print "Reading in fort.q**** files"

    ## Read in fort.q**** files
    outdir = 'run%s/_output/' %str(run).zfill(2)
    for fnum in fnums:
        fname = outdir + 'fort.q%s' % str(fnum).zfill(4)
        varname = 'fort_%s' % str(fnum).zfill(4)
#        print 'Reading in: '+ fname
#        print 'varname is: '+ varname
        exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))
        status(fnum,nfiles)

        exec('y = %s[:,%r].copy()'%(varname,measure_num))
        z = y.reshape(n,m,order='F')
        # Corrects for the fact h is in m and eta is in cm
        if measure_num == 3:
            z = z
        else:
            z = z

        z = z*1./100
        if fnum == 1:
#            print "saving z1 at k = %s" %k
            z1 = z
        max_h = where(z>max_h,z,max_h)  #computes max over entire domain
    return(max_h,z1)                            

################
## Status Bar ##
################

def status(n,N):
    n = float(n)
    N = float(N)
    percent = n/N*100
    sys.stdout.write("[==> ]%3d%%\r" %percent)
    sys.stdout.flush()


####################
## mu Calculation ##
####################
def muCalc(eta,zeta_i,fieldType,nx,ny,nu,T):
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
    
    for s in range(0,runs):         #loop through events
        for ix in range(0,nx):     #loop through x
            status(ix,nx)
            for iy in range(0,ny): #loop through y
                #The following equations are from Gonzalez et al 2009
                zeta_0j = eta[ix,iy] + MSL + C*(MHHW-MSL)*exp(-alpha*(eta[ix,iy]/sigma_0)**beta) #eqn 2b
                sigma_j = sigma_0-CP*sigma_0*exp(-1*alphaP*(eta[ix,iy]/sigma_0)**betaP)          #eqn 2c
                mu[ix,iy] = mu[ix,iy] + 1./2*nu[s]*(1-erf((zeta_i-zeta_0j)/(sqrt(2)*sigma_j)))   #eqn 3b/5

  
    return(mu)

###################
## P Calculation ##
###################

def pCalc(mu,nx,ny,T):
    print "Computing P(zeta>zeta_i)"
    P = zeros((nx,ny))
    P = ones((nx,ny))-exp(-1*mu*T) #eqn 6
    return(P)
####################
## Plotting Tools ##
####################
def plotting():

    #masking nessisary arrays
    mask_lim = -10
    masked_max_h = numpy.ma.masked_where(z1 > mask_lim,maxNum)
    masked_P = numpy.ma.masked_where(z1 > mask_lim,P)
    masked_mu = numpy.ma.masked_where(z1 > mask_lim,mu)
    coast = numpy.ma.masked_where(z1 > 0, ones((nx,ny)))
    h_non_zero = numpy.ma.masked_where(maxNum == 0,maxNum)
    zInit = numpy.ma.masked_where(z1 < -10**10, z1)
    #calculating axis corrections:
    dx = 2.777778*10**(-4)   #latitude dx
    dy = 2.777778*10**(-4)  #longitude dy
    x0 = 235.7655#*10**2      #lower left latitude of the grid
    y0 = 41.76956#*10**1      #lower left longitude of the grid
    x = linspace(x0,x0+dx*nx,nx)
    y = linspace(y0+dy*ny,y0,ny)

    figure(1)
    clf()
    pcolormesh(x,y,rot90(masked_max_h))
    axis([x[0], x[-1], y[-1], y[0]])
    jet()
    colorbar()
    title('$P_i(\zeta>\zeta_i)$')
    savefig('%sProbabilityPlot.png' %measure[measure_num])
    
    figure(2)
    clf()
    pcolormesh(x,y,rot90(masked_mu))
    axis([x[0], x[-1], y[-1], y[0]])
    jet()
    colorbar()
    title('$\mu_{i,j}$')
    savefig('%sMaxMaskedMuPlot.png'%measure[measure_num])

    figure(3)
    clf()
    pcolormesh(x,y,rot90(masked_max_h))
    axis([x[0], x[-1], y[-1], y[0]])
    jet()
    colorbar()
    title('$\hat{\eta}_i$ w.r.t. Topography')
    savefig('%sMaxPlotMasked.png' %measure[measure_num] )
    
    figure(4)
    clf()
    pcolormesh(x,y,rot90(coast))
    axis([x[0], x[-1], y[-1], y[0]])
    summer()
    #   colorbar()
    title('Crescent City Coast')
    savefig('%sCoastLine.png' %measure[measure_num])

    figure(5)
    pcolormesh(x,y,rot90(h_non_zero))
    axis([x[0], x[-1], y[-1], y[0]])
    jet()
    colorbar()
    title('Where does flooding occur')
    savefig('%sFloodWhere.png' %measure[measure_num])
    
    figure(6)
    clf()
    pcolormesh(x,y,rot90(zInit))
    axis([x[0], x[-1], y[-1], y[0]])
    summer()
    colorbar()
    title('Crescent City Coast')


    show()
    return(coast)


##################
## Main Program ##
##################


#test = 0 implies to run code on full set of  inundations
#test = 1 implies to run a test on a 2x2 identity matrix instead of inundation

test = 0
#test case with simple runup cases
if test == 1:
    zeta_i = 2.1                       #meters of inundation
    fieldType = 1                      #1 = far-field, 2 = near-field
    ((nx,ny,runs)) = ((2,2,1))
    field_1 = array([[1,0],[0,1]])*100
    field_2 = array([[1,1],[0,0]])*100

    T_M = 520
    T = 1
    nu = ones((1,runs))*1./T_M

    (P,mu_1) = tidalUncert(field_1,zeta_i,fieldType,nx,ny,runs,nu,T)
    (P,mu_2) = tidalUncert(field_2,zeta_i,fieldType,nx,ny,runs,nu,T)

    mu = mu_1 + mu_2 
    P = ones((nx,ny))-exp(-1*mu*T)
#use real data
if test == 0:
    measure = (('H','HU','HV','ETA'))  #the parameter to measure
    measure_num = 0                    #identifying number of measure
    zeta_i = 2.1                       #meters of inundation
    fieldType = 1                      #1 = far-field, 2 = near-field
    T_M = 520                          #Recurance Time
    T = 1                              #Time interval of interest

    zeta_i = 2.1                       #meters of inundation
    fieldType = 1                      #1 = far-field, 2 = near-field
    nfiles = 59
    ((nx,ny,runs)) = ((271,192,1))
    nu = ones((1,runs))*1./T_M
    mu = zeros((nx,ny))
    muRun = zeros((nx,ny))
    runList = range(15,18)
    for run in runList:
        maxNum = 'max_h_%s' % str(run).zfill(2)
        z1Num = 'z1_%s' % str(run).zfill(2)
        (maxNum,z1Num) = maxH(nfiles,run,nx,ny,runs)
        muRun = muCalc(maxNum,zeta_i,fieldType,nx,ny,nu,T) 
        mu = mu + muRun
    P = pCalc(mu,nx,ny,T)
    z1 = z1Num
    plotting()
