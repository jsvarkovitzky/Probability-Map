#his program reads in the fort.q**** files and generates a max inundation 
#Jonathan Varkovitzky


from numpy import *
from matplotlib import *
from matplotlib import rc
rc("text",usetex=True)
from pylab import *
from scipy.special import *
#inundationMax()

########################################################
## Iterate over each time step to find max inundation ##
########################################################
def maxH( nfiles,n,m,k):
    fnums = range(1,nfiles)
    max_h = zeros((n,m))
    print "Reading in fort.q**** files"

    ## Read in fort.q**** files
    for fnum in fnums:
        outdir = 'MOST_TEST/_output/'
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
            print "saving z1 at k = %s" %k
            z1 = z
        max_h = where(z>max_h,z,max_h)  #computes max over entire domain
    return(max_h,z1)                            

################################################
## Take max_h information and generate P_{ij} ##
################################################
def tidalUncert(eta,zeta_i,fieldType,nx,ny,runs,nu,T):
    print "Computing P(zeta>zeta_i)"

    
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
    
    #calculate the probability of excedence at each point (i,j)
    mu = zeros((nx,ny))
    P = zeros((nx,ny))
    
    print 'Computing Probabilities: '
    for s in range(0,runs):         #loop through events
        for ix in range(0,nx):     #loop through x
            status(ix,nx)
            for iy in range(0,ny): #loop through y
                #The following equations are from Gonzalez et al 2009
                zeta_0j = eta[ix,iy] + MSL + C*(MHHW-MSL)*exp(-alpha*(eta[ix,iy]/sigma_0)**beta) #eqn 2b
                sigma_j = sigma_0-CP*sigma_0*exp(-1*alphaP*(eta[ix,iy]/sigma_0)**betaP)          #eqn 2c
                mu[ix,iy] = mu[ix,iy] + 1./2*nu[s]*(1-erf((zeta_i-zeta_0j)/(sqrt(2)*sigma_j)))   #eqn 3b/5

    P = ones((nx,ny))-exp(-1*mu*T)                                                               #eqn 6
    return(P,mu)

####################
## Plotting Tools ##
####################
def plotting():

    #masking nessisary arrays
    mask_lim = -10
    masked_max_h = numpy.ma.masked_where(z1 > mask_lim,max_h)
    masked_P = numpy.ma.masked_where(z1 > mask_lim,P)
    masked_mu = numpy.ma.masked_where(z1 > mask_lim,mu)
    coast = numpy.ma.masked_where(z1 > 0, ones((nx,ny)))
    h_non_zero = numpy.ma.masked_where(max_h == 0,max_h)
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
## Program Main ##
##################
"""
close ('all')
#test = 0 implies to run code on full set of  inundations
#test = 1 implies to run a test on a 2x2 identity matrix instead of inundation
test = 0
measure = (('H','HU','HV','ETA'))  #the parameter to measure
measure_num = 0                    #identifying number of measure
zeta_i = 2.1                       #meters of inundation
fieldType = 1                      #1 = far-field, 2 = near-field


## Testing max_h generation ##

if test == 0:
    nfiles = 59
    ((nx,ny,runs)) = ((271,192,1)) #x-pts, y-pts, k-timesteps ### in the future automate this!!!!!
    #compute raw max height from simulation
    (max_h,z1) = maxH(nfiles,nx,ny,runs)
  

## Actual Run to calculate max_h ##
  
if test == 1:
    ((nx,ny,runs)) = ((2,2,1))
    max_h = array([[1,0],[0,1]])*100
    z1 = max_h



T_M = 520                 #recurance time from Gonzalez et al. 2009
T = 1                   #period of intrest 
nu = ones((1,runs))*1./T_M


(P,mu) = tidalUncert(max_h,zeta_i,fieldType,nx,ny,runs,nu,T)

#Calls plotting algorithms
plotting()
"""
