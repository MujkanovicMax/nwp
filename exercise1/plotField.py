#!/usr/bin/env python
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap


#Read reanalysis data
filename='eraint_2019010100.nc'
ncf = netCDF4.Dataset(filename, 'r')
#prints the content
print(ncf)
lons=ncf.variables['longitude'][:]
lats=ncf.variables['latitude'][:]
u=ncf.variables['u'][0,:]
v=ncf.variables["v"][0,:]
vo=ncf.variables["vo"][0,:]
time=ncf.variables["time"][0]
z=ncf.variables["z"][0,:]

varlist = ['u',"v","vo","z"]

ind_x = np.where(lats>=0)
lats = lats[ind_x[0]]

coord = np.meshgrid(lons,lats)



for i,var in enumerate(varlist):
    
    map = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
    map.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
    map.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
    map.contourf(coord[0],coord[1],ncf[var][0,:][ind_x[0]])
    
    plt.tight_layout()
    
    plt.savefig("./field3_" + var +".png")
    
    
ncf.close()   

#ind_x = np.where(lats > 30 and lats < 60)


#mean_zonal_wind = np.nanmean(u)

#radius_e = 6371000
#dx = radius_e*np.pi/180.* np.abs(np.cos(lats))
#dy = radius_e*np.pi/180.

#w2 = vo*vo
#enstrophy = np.sum(np.sum(w2*dy,axis=1)*dx,axis=0)/(4*np.pi*radius_e*radius_e)



