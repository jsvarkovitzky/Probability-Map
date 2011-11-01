#This program reads in the fort.q**** files and generates a max inundation 
#Jonathan Varkovitzky


from numpy import *
from matplotlib import *
from pylab import *

#inundationMax()

########################################################
## Iterate over each time step to find max inundation ##
########################################################
def maxH( nfiles,n,m ):
    fnums = range(1,nfiles)
    print "Reading in fort.q**** files"

    ## Read in fort.q**** files
    for fnum in fnums:
        outdir = 'MOST_TEST/_output/'
        fname = outdir + 'fort.q%s' % str(fnum).zfill(4)
        varname = 'fort_%s' %str(fnum).zfill(4)
        #print 'Reading in: '+ fname
        exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))
    
        max_h = zeros((n,m))
        max_h2 = zeros((n,m))
    
    print "Analysing frames"
    for k in range(1,nfiles):
        varnam = 'fort_%s' %str(k).zfill(4)
        exec('y = %s[:,1].copy()'%varname)
        z = y.reshape(n,m,order='F')
        z = z*1./100
        if k == 1:
            z1 = z
            
        #print 'Analysing frame # %s' %str(k).zfill(4)
        for i in range(1,n):
            for j in range(1,m):
                if  z[i,j] < 10**20 and z[i,j]>max_h[i,j]:
                     max_h[i,j] = z[i,j]
                        
    return(max_h)                            

################################################
## Take max_h information and generate P_{ij} ##
################################################
def tidalUncert(zeta,zeta_i,fieldType,n,m):
    print "Computing P(zeta>zeta_i)"
    #Constant values obtained from Mofjeld et al. 2007
    MSL = 1.130
    MHHW = 2.095
    alpha = 0.170
    beta = 0.858
    C = 1.044
    CPrime = 0.707 
    alphaPrime = 0.056
    betaPrime = 1.119
    sigma_0 = 0.638 
    
    #calculate the probability of excedence at each point (i,j)
    for i in range(1,n):
        for j in range(1,m): 
            zeta_0 = zeta[i,j] + MSL + C*(MHHW-MSL)*exp(-alpha*(zeta[i,j]/sigma_0)**betaPrime)

    return


##################
## Program Main ##
##################

nfiles = 6
((n,m)) = ((271,192)) #x-pts, y-pts ### in the future automate this!!!!!
#compute raw max height from simulation
max_h = maxH(nfiles,n,m)
#compute P_{ij}
zeta_i = 1              #meters of inundation
fieldType = 1           #1 = far-field, 2 = near-field
P = tidalUncert(max_h,zeta_i,fieldType,n,m)

#########################
## Plotting Algorithms ##
#########################

figure(3)
pcolormesh(max_h)
colorbar()
show()





#plot_where = numpy.ma.masked_where(z1>10**20,max_h)

#figure(1)
#matplotlib.pyplot.contour(max_h,100)
#colorbar()

#figure(2)
#matplotlib.pyplot.contour(plot_where,100)
#colorbar()




