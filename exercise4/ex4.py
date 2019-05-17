#!/usr/bin/env python
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap



def ddx(values,dx):
    dvalues = np.zeros(values.shape)
    dvalues[0,:]      = values[-1,:]-values[1,:] 
    dvalues[1:-1,:]   = values[0:-2,:] - values[2:,:]
    dvalues[-1,:]     = values[-2,:] - values[0,:]
    
    dvalues = dvalues / (2*dx)
    
    return dvalues

def ddy(values,dx):
    dvalues = np.zeros(values.shape)
    dvalues[:,0]      = (values[:,0]-values[:,1])/dx 
    dvalues[:,1:-1]   = (values[:,0:-2] - values[:,2:])/(2*dx)
    dvalues[:,-1]     = (values[:,-2] - values[:,-1])/dx
    
    return dvalues


#Read reanalysis data
filename='eraint_2019020100.nc'
ncf = netCDF4.Dataset(filename, 'r')
filename='eraint_2019020106.nc'
ncf2 = netCDF4.Dataset(filename, 'r')
u_an=ncf2.variables['u'][0,:]
v_an=ncf2.variables["v"][0,:]

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
u = u[ind_x]
v = v[ind_x]
v_an=v_an[ind_x]
u_an=u_an[ind_x]


lats=lats[::-1]
u = u[::-1,:]
v = v[::-1,:]
u_an = u_an[::-1,:]
v_an = v_an[::-1,:]

coord = np.meshgrid(lons,lats)
radius_e = 6371000
dx = radius_e*np.pi/180.* np.cos(45/180.*np.pi)
dy = radius_e*np.pi/180.
dt = 6*3600.

vorticity = np.zeros(v.shape)
tmp = vorticity*1
vorticity_an = vorticity*1
tmp2 = tmp*1

#for i in range(v.shape[1]):
    #vorticity[:,i] = ddx(v[:,i],dx)
    #vorticity_an[:,i] = ddx(v_an[:,i],dx)
    
#for j in range(v.shape[0]):
     #tmp[j,:] = - ddy(u[j,:],dy)
     #tmp2[j,:] = -ddy(u_an[j,:],dy)
vorticity=ddx(v,dx)-ddy(u,dy)     
     

     
vorticity = vorticity + tmp 
vorticity_an = vorticity_an + tmp2

#for i in range(vorticity.shape[1]):
    #tmp[:,i] = ddx(vorticity[:,i],dx)
    
#for j in range(vorticity.shape[0]):
    #tmp2[j,:] = - ddy(vorticity[j,:],dy)

vorticity2 = vorticity - ((u*tmp + v*tmp2) + v*1.6e-11)*dt

    

map = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map.contourf(coord[0]-180.,coord[1],vorticity_an[ind_x[0]])
    
plt.tight_layout()
    
plt.savefig("test_an.png")

map2 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map2.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map2.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map2.contourf(coord[0]-180.,coord[1],vorticity2[ind_x[0]])
    
plt.tight_layout()
    
plt.savefig("test_fc.png")


map2 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map2.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map2.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map2.contourf(coord[0]-180.,coord[1],vorticity[ind_x[0]])
    
plt.tight_layout()
    
plt.savefig("test_init.png")