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

             
#    print "Analysing frames"
#    print "Analysing Frames: "
#    for k in range(0,nfiles):
#        status(k,nfiles)
#        varname = 'fort_%s' % str(k).zfill(4)
#        print varnam
        exec('y = %s[:,3].copy()'%varname)
        z = y.reshape(n,m,order='F')
        z = z*1./100
        
        if fnum == 1:
            print "saving z1 at k = %s" %k
            z1 = z

        if fnum == 10:
            z2 = z
            #print z2 == z1

        #check and correct max_h against frame k
        for i in range(0,n):
            for j in range(0,m):
                if z[i,j]>max_h[i,j]:
                     max_h[i,j] = z[i,j]
#    figure(4)
#    pcolormesh(z1)
#    colorbar()
#    show()
                        
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
    
    print 'Computing Probability: '
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
    
    figure(1)
    pcolormesh(rot90(P))
    colorbar()
    title('Plot of P')
    savefig('ProbabilityPlot.png')

    figure(2)
    pcolormesh(rot90(max_h))
    colorbar()
    title('Plot of max_h')
    savefig('MaxHeightPlot.png')

    #Plotting Masked Plots:

    masked_max_h = numpy.ma.masked_where(z1 <= 0,max_h)
    coast = numpy.ma.masked_where(z1 < 0, ones((nx,ny)), )
    
    figure(3)
    pcolormesh(rot90(masked_max_h))
    colorbar()
    title('masked plot of max_h')
    savefig('MaxHeightPlot.png')

    figure(4)
    pcolormesh(rot90(coast))
    colorbar()

    figure(5)
    pcolormesh(z1)
    colorbar()
    show()


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

#test = 0 implies to run code on full set of  inundations
#test = 1 implies to run a test on a 2x2 identity matrix instead of inundation
test = 0

if test == 0:
    nfiles = 59
    ((nx,ny,runs)) = ((271,192,1)) #x-pts, y-pts, k-timesteps ### in the future automate this!!!!!
    #compute raw max height from simulation
    (max_h,z1) = maxH(nfiles,nx,ny,runs)
    
if test == 1:
    ((nx,ny,runs)) = ((2,2,1))
    max_h = array([[1,0],[0,1]])*100
    z1 = max_h

T_M = 520                 #recurance time from Gonzalez et al. 2009
T = 1                   #period of intrest 
nu = ones((1,runs))*1./T_M

#compute P_{ij}
zeta_i = 2.1            #meters of inundation
fieldType = 1           #1 = far-field, 2 = near-field
(P,mu) = tidalUncert(max_h,zeta_i,fieldType,nx,ny,runs,nu,T)

#Calls plotting algorithms
coast = plotting()






#plot_where = numpy.ma.masked_where(z1>10**20,max_h)

#figure(1)
#matplotlib.pyplot.contour(max_h,100)
#colorbar()

#figure(2)
#matplotlib.pyplot.contour(plot_where,100)
#colorbar()




