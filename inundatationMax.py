#This program reads in the fort.q**** files and generates a max inundation 
#Jonathan Varkovitzky
from numpy import *



##############################
## Read in fort.q**** files ##
##############################
nfiles = 5
fnums = range(1,nfiles)

for fnum in fnums:
    outdir = 'MOST_TEST/_output/'
    fname = outdir + 'fort.q%s' % str(fnum).zfill(4)
    varname = 'fort_%s' %str(fnum).zfill(4)
    print 'Reading in: '+ fname
    exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))

########################################################
## Iterate over each time step to find max inundation ##
########################################################
((n,m)) = ((271,192)) #x-pts, y-pts ### in the future automate this!!!!!
max_h = zeros((n,m))

for k in range(1,nfiles):
    varnam = 'fort_%s' %str(k).zfill(4)
    exec('y = %s[:,1].copy()'%varname)
    z = y.reshape(n,m,order='F')
    for i in range(1,n):
        for j in range(1,m):
            if z[i,j]>max_h[i,j]:
                max_h[i,j] = z[i,j]



