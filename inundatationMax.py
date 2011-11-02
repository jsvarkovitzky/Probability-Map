#This program reads in the fort.q**** files and generates a max inundation 
#Jonathan Varkovitzky


from numpy import *
from matplotlib import *
from pylab import *
from scipy.special import *
#inundationMax()

########################################################
## Iterate over each time step to find max inundation ##
########################################################
def maxH( nfiles,n,m,k ):
    fnums = range(1,nfiles)
    print "Reading in fort.q**** files"

    ## Read in fort.q**** files
    for fnum in fnums:
        outdir = 'MOST_TEST/_output/'
        fname = outdir + 'fort.q%s' % str(fnum).zfill(4)
        varname = 'fort_%s' %str(fnum).zfill(4)
        #print 'Reading in: '+ fname
        exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))
        status(fnum,nfiles)
        max_h = zeros((n,m))
        max_h2 = zeros((n,m))
    
    print "Analysing frames"
    for k in range(0,nfiles):
        status(k,nfiles)
        varnam = 'fort_%s' %str(k).zfill(4)
        exec('y = %s[:,1].copy()'%varname)
        z = y.reshape(n,m,order='F')
        z = z*1./100
        if k == 1:
            z1 = z
            
        #print 'Analysing frame # %s' %str(k).zfill(4)
        for i in range(0,n):
            for j in range(0,m):
                if  z[i,j] < 10**20 and z[i,j]>max_h[i,j]:
                     max_h[i,j] = z[i,j]
                        
    return(max_h)                            

################################################
## Take max_h information and generate P_{ij} ##
################################################
def tidalUncert(zeta,zeta_i,fieldType,n,m,k,nu,T):
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
    mu = zeros((n,m))
    P = zeros((n,m))
    #diff = zeros((n,m))
    for s in range(0,k):         #loop through events
        for i in range(0,n):     #loop through x
            for j in range(0,m): #loop through y
                #equations from Gonzalez et al 2009
                zeta_0 = zeta[i,j] + MSL + C*(MHHW-MSL)*exp(-alpha*(zeta[i,j]/sigma_0)**beta)
                sigma = sigma_0-CP*sigma_0*exp(-1*alphaP*(zeta[i,j]/sigma_0)**betaP)
                mu[i,j] = 1./2*nu[s]*(1-erf((zeta_i-zeta[i,j])/(sqrt(2)*sigma)))
                print "(i,j) = (%s,%r)" %(i,j)
                print "The vale of zeta_0 is %s" %(zeta_0)
                print "The value of sigma is %s" %(sigma)
                print "The value of mu[i,j] is %s" %(mu[i,j])
    P = ones((n,m))-exp(-1*mu*T)
    return(P,mu,diff)

####################
## Plotting Tools ##
####################
def plotting():
    figure(1)
    pcolormesh(P)
    colorbar()
    show()

################
## Status Bar ##
################

def status(n,N):
    n = float(n)
    N = float(N)
    percent = n/N*100
    sys.stdout.write("[==> ]%3d%%\r" % percent)
    sys.stdout.flush()

##################
## Program Main ##
##################

nfiles = 59
((n,m,k)) = ((271,192,1)) #x-pts, y-pts, k-timesteps ### in the future automate this!!!!!
#((n,m,k)) = ((2,2,1))
T_M = 520                 #recurance time from Gonzalez et al. 2009
T = 100                   #period of intrest 
nu = ones((1,k))*1./T_M

#compute raw max height from simulation
max_h = maxH(nfiles,n,m,k)

#max_h = array([[1,0],[0,1]])
#compute P_{ij}
zeta_i = 0.0              #meters of inundation
fieldType = 1           #1 = far-field, 2 = near-field
(P,mu,diff) = tidalUncert(max_h,zeta_i,fieldType,n,m,k,nu,T)

#Calls plotting algorithms
plotting()






#plot_where = numpy.ma.masked_where(z1>10**20,max_h)

#figure(1)
#matplotlib.pyplot.contour(max_h,100)
#colorbar()

#figure(2)
#matplotlib.pyplot.contour(plot_where,100)
#colorbar()




