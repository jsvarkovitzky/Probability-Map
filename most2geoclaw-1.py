"""
Module for converting MOST data and output to GeoClaw format.
"""

import os, glob, re, sys, numpy
from scipy.io.netcdf import netcdf_file
import matplotlib.pyplot as plt
from gauge_info import gauge_info as gauge_info
from numpy import ones as ones
from matplotlib.mlab import griddata
from numpy import *
from scipy import interpolate
class Most2Geoclaw(object):

#    def __init__(self,fname = 'CC'):    
    def __init__(self,fname = 'source15'):
        h = netcdf_file(fname + '_ha.nc','r')
        u = netcdf_file(fname + '_ua.nc','r')
        v = netcdf_file(fname + '_va.nc','r')


#        h = netcdf_file(fname + '_runup_ha.nc','r')
#        u = netcdf_file(fname + '_runup_ua.nc','r')
#        v = netcdf_file(fname + '_runup_va.nc','r')

        self.ha = h.variables['HA'][:]
        self.ua = u.variables['UA'][:]
        self.va = v.variables['VA'][:]

        self.lat = h.variables['LAT'][:]
        self.lat = list(self.lat)
        self.backwards = True
        if self.backwards == True:
            self.lat.reverse()
        self.funky_frank = True
#        print 'len self.lat = ', len(self.lat)
        self.lon = h.variables['LON'][:]
#        print 'len self.on = ',len(self.lon)
        self.xll = self.lon[0]
        self.dx = self.lon[1] - self.xll
        self.yll = self.lat[0]
        self.dy = self.lat[1] - self.yll
        self.gauge_no = len(gauge_info)

        self.nrows = self.ha[0].shape[0]
        self.ncols = self.ha[0].shape[1]
        self.deltat=1

        # outputs are given as time - row - column nested arrays

        self.time = h.variables['TIME'][:]
        
        self.gauge_lon = []
        self.gauge_lat = []
        for gauge in range(len(gauge_info)):
            self.gauge_lon.append(gauge_info[gauge][0])
            self.gauge_lat.append(gauge_info[gauge][1])


    def most2tt3(self,fname):
    
    #    Converts MOST topo file to tt3 format.

        f = open(fname).readlines()
        mn = f[0].split()
        ncols = int(mn[0])
        nrows = int(mn[1])
        xll = float(f[1]) 
        dx = float(f[2]) - xll
        xll = xll - 360.
        yll = float(f[nrows+ncols])
        dy =  float(f[nrows+ncols-1]) - yll
#        if abs(dx-dy) > 1.e-6:
#            print '*** WARNING: dx = ',dx,'  dy = ',dy
        cellsize = dx

        return f[nrows+ncols+1:]






    def most2fortt(self,fname,t):

    #        Converts MOST output files to fort.t files.

            fortname = 'fort.t' + str(t).zfill(4)
            f = open('./_output/'+fortname,'w')

            f.write("%18.8e     time\n" % self.time[t])
            f.write("%5i                  meqn\n" % 4)
            f.write("%5i                  ngrids\n" % 1)
            f.write("%5i                  naux\n" % 3)
            f.write("%5i                  ndim\n" % 2)

            f.close()
            print "Created %s from %s at time t = %s" % (fortname, fname, t)


    def most2fortq(self,fname,t):

    #        Converts MOST output files to fort.q files.

            
    #        if abs(dx-dy) > 1.e-6:
    #            print '*** WARNING: dx = ',dx,'  dy = ',dy
    #        cellsize = dx

            fortname = 'fort.q' + str(t).zfill(4)
            f2 = open('./_output/'+fortname,'w')
            f = open(bathyname).readlines()
            f2.write("%5i                  grid_number\n" % 1)
            f2.write("%5i                  AMR_level\n" % 1)
            f2.write("%5i                  mx\n" % self.ncols)
            f2.write("%5i                  my\n" % self.nrows)
            f2.write("%8e                  xlow\n" % self.xll)
            f2.write("%8e                  ylow\n" % self.yll)
            f2.write("%8e                  dx\n" % self.dx)
            f2.write("%8e                  dy\n" % self.dy)
            f2.write("\n")

            for r in range(len(self.ha[t])):
                if self.backwards == True: 
                    increment = -1-r            #this is because clawpack and most save model data differently
                else: 
                    increment = r
                k = p.most2tt3(bathyname)
                bathy = k[(increment)].split()

                if self.funky_frank == True:
                    increment = -1-increment   # indexes the eta values backwards
#                plotter(ha[t])
                for c in range(len(self.ha[t][increment])):
                    bathy_value = -float(bathy[c])
                    eta = self.ha[t][increment][c]/100
                    h = eta - bathy_value
                    hu = self.ua[t][increment][c]/100*h
                    hv = self.va[t][increment][c]/100*h

                    f2.write("%18.8e    " % h)
                    f2.write("%18.8e    " % hu)
                    f2.write("%18.8e    " % hv)
                    f2.write("%18.8e    " % eta)
                    f2.write("\n")
                f2.write("\n")
            f2.close()
            print "Created %s from %s" % (fortname,fname)
        
            yul = self.yll + self.nrows  * self.dy
            xul = self.xll + self.ncols * self.dx
            f3 = open('./plotvar.py','w')
            f3.write("plotvar = {'xll': %8e, 'yll': %8e, 'xul':%8e, 'yul':%8e}" %(self.xll,self.yll,xul,yul))
        
    def most2fortgauge(self,fname,t,gauge_no):

#        Converts MOST output files to fort.q files.


#        if abs(dx-dy) > 1.e-6:
#            print '*** WARNING: dx = ',dx,'  dy = ',dy
#        cellsize = dx

        fortname = 'fort.gauge'
        f2 = open('./_output/'+fortname,'a')
        f = open(bathyname).readlines()
        gauge_time = self.time[0] + t*self.deltat
        bathy = p.most2tt3(bathyname)
        bathy_array= p.bathy_to_array(bathy)

        for gauge in range(len(gauge_info)):
            gauge_id = gauge_info[gauge][2]
            xcoord = p.find_grid_points(self.gauge_lon[gauge],self.xll,self.dx)
            ycoord = p.find_grid_points(self.gauge_lat[gauge],self.yll,self.dy)
#            print gauge
#            print 'self.gauge_lon[gauge] = ', self.gauge_lon[gauge]
#            print 'self.gauge_lat[gauge] = ', self.gauge_lat[gauge]
#            print 'ycoord = ', ycoord
#            print 'xcoord = ', xcoord
            ycoord = -1 - ycoord
            eta = self.ha[t][ycoord][xcoord]/100
            h = eta + bathy_array[ycoord][xcoord]          # eta - (-bathy)  because bathy is defined in other direction for clawpack
            hu = self.ua[t][ycoord][xcoord]/100*h
            hv = self.va[t][ycoord][xcoord]/100*h

            
            f2.write("   %6.0i  " % gauge_id)    #gauge id number
            f2.write("%1.0i " % 1)           #grid number, will always be 1 for most models
            f2.write("%18.8e " % gauge_time) #print out the gauge information at every time step
            f2.write("%18.8e " % h)
            f2.write("%18.8e " % hu)
            f2.write("%18.8e " % hv)
            f2.write("%18.8e" % eta)
            f2.write("\n")
        f2.close()
        print "Created %s from %s" % (fortname,fname)

    def bathy_to_array(self,bathy):       #   feeds in a list of strings - outputs an array
        nrows_bathy = len(bathy)
        ncols_bathy = len(bathy[0].split())
        bathy_array = ones((nrows_bathy,ncols_bathy))
        for k in range(len(bathy)):
            bathy_list_of_strings = bathy[k].split()
            for d in range(len(bathy_list_of_strings)):
                bathy_array[k][d] = bathy_array[k][d]*float(bathy_list_of_strings[d])
        return bathy_array
            
    def find_grid_points(self,coord,axisll,delta):
        scaledcoor = coord - axisll
        gridpoint = round(scaledcoor/delta)
        gridpoint = int(gridpoint)
        return gridpoint
    
    
    def most2setgauges(self,gauge_no,dir):
        fortname = 'setgauges.data'
        f2 = open(dir+fortname,'a')
#        f2 = open('./'+fortname,'a')
        f2.write("########################################################") 
        f2.write("\n")
        f2.write("### DO NOT EDIT THIS FILE:  GENERATED AUTOMATICALLY ####") 
        f2.write("\n")
        f2.write("### To modify data, edit  most2geoclaw.py           ####") 
        f2.write("\n")
        f2.write("########################################################") 
        f2.write("\n")
        f2.write("########################################################") 
        f2.write("\n")
        f2.write("\n")
        f2.write("%4i   =: ngauges" % gauge_no)
        f2.write("\n")
        for gauge in range(gauge_no):
            gaugename = int(gauge)+1
            f2.write("%4i" % gaugename)
            f2.write("  %18.8e" % self.gauge_lon[gauge])
            f2.write("   %18.8e" % self.gauge_lat[gauge])
            f2.write("   %18.8e" % self.time[0])     
            f2.write("   %18.8e" % self.time[-1])
            f2.write("\n")
        f2.close
        print "Created %s" % (fortname)
#    def plotter(self,data):
#            
#            x , y = meshgrid(self.lat, self.lon)
#            m = Basemap()
#            m.contourf(x, y, data)
#            plt.show()
        
    def builder(self,fname):
        f2 = open('./_output/fort.gauge','w')
        f2.close
        f2 = open('./setgauges.data','w')
        f2.close
        f2 = open('./_output/setgauges.data','w')
        f2.close
        for t in range(0,len(self.time),self.deltat):    
            p.most2fortq(fname,t)
            p.most2fortt(fname,t)
            p.most2fortgauge(fname,t,self.gauge_no)
        dir = './'
        p.most2setgauges(self.gauge_no,dir)
        dir = './_output/'
        p.most2setgauges(self.gauge_no,dir)  
#        print 'self.xll = ', self.xll 
#        print 'self.yll = ', self.yll
#        print 'self.lat[0] = ', self.lat[0]
#        print 'self.lat[-1] = ', self.lat[-1]
#        print 'self.lon[0] = ', self.lon[0]
#        print 'self.lon[-1] = ', self.lon[-1]          
#        print 'self.dy = ', self.dy
#        print 'self.dx = ', self.dx
            
if __name__=='__main__':
    import sys
    import os
#    fname = sys.argv[1]
#    fname = 'CC'
#    bathyname = 'cresc1secm_mod.txt'
    fname = 'source15'
    bathyname = 'cresc1secm_mod.asc.s.c'
    p = Most2Geoclaw()
    p.builder(fname)
#    p.plotter()
#    most2tt3(bathyname)
    
