#This program tests reading in n files into various arrays
#Jonathan Varkovitzky
from numpy import *

fnums = range(1,60)

for fnum in fnums:
    outdir = 'MOST_TEST/_output/'
    fname = outdir + 'fort.q%s' % str(fnum).zfill(4)
    varname = 'fort_%s' %str(fnum).zfill(4)
    print fname
    exec('%s=loadtxt(%r, skiprows=9)'%(varname,fname))
