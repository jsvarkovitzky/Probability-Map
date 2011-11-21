"""
Module for converting GeoClaw files into MOST format.
"""

import os, glob, re
from numpy import *

def tt32most(fname):
    print fname
    f = open(fname).readlines()
    ncols = int(f[0].split()[0])
    nrows = int(f[1].split()[0])
    xll = float(f[2].split()[0])+360
    yll = float(f[3].split()[0])
    dx = float(f[4].split()[0])

#    print ncols, nrows, xll, yll, dx
    topoData = loadtxt(fname, skiprows=6)
    topoData = rot90(transpose(topoData))
    fname2 = os.path.splitext(fname)[0] + '.2c'
    f2 = open(fname2,'w')
    f2.write(' %s %r\n' %(ncols, nrows))

    for i in range(0,ncols):
        x = xll+dx*i
        f2.write('   %s\n' %x )
    for j in range(0,nrows):
        y = yll+dx*(nrows - (j+1))
        f2.write('   %s\n' %y )
    for k in range(6,nrows):
        f2.write('%s\n' %f[k])

    f2.close()
    print 'Created ',fname2

#first write out 271 y values then 192 x values, then topo data
    
if __name__ =='__main__':
    import sys
    tt32most(sys.argv[1])

