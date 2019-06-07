#!/usr/bin/env python
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats


def ddx(values,dx):
    dvalues = np.zeros(values.shape)
    dvalues[0,:]      = -values[-1,:]+values[1,:] 
    dvalues[1:-1,:]   = -values[0:-2,:] + values[2:,:]
    dvalues[-1,:]     = -values[-2,:] + values[0,:]
    
    dvalues = dvalues / (2*dx)
    
    return dvalues

def ddy(values,dx):
    dvalues = np.zeros(values.shape)
    dvalues[:,0]      = (-values[:,0]+values[:,1])/dx 
    dvalues[:,1:-1]   = (-values[:,0:-2] + values[:,2:])/(2*dx)
    dvalues[:,-1]     = (-values[:,-2] + values[:,-1])/dx
    
    return dvalues


#Read reanalysis data
filename='../exercise4/eraint_2019020100.nc'
ncf = netCDF4.Dataset(filename, 'r')
filename='../exercise4/eraint_2019020106.nc'
ncf2 = netCDF4.Dataset(filename, 'r')
u_an=ncf2.variables['u'][0,:]
v_an=ncf2.variables["v"][0,:]
vo_an=ncf2.variables["vo"][0,:]

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
vo = vo[ind_x]
vo_an = vo_an[ind_x]


lats=lats[::-1]
u = u[::-1,:].transpose()
v = v[::-1,:].transpose()
u_an = u_an[::-1,:].transpose()
v_an = v_an[::-1,:].transpose()
vo = vo[::-1,:].transpose()
vo_an = vo_an[::-1,:].transpose()

coord = np.meshgrid(lons-180.,lats)
radius_e = 6371000
dx = radius_e*np.pi/180.* np.cos(45/180.*np.pi)
dy = radius_e*np.pi/180.
dt = 6*3600.

#vorticity = np.zeros(v.shape)
#vorticity=ddx(v,dx)-ddy(u,dy)     


#vorticity2 = vorticity - (u*ddx(vorticity,dx) + v*ddy(vorticity,dy) + v*1.6e-11)*dt

vo_en = np.hstack((np.fliplr(vo[:,1:-1])*(-1),vo))
vo_en[:,0] = 0
vo_en[:,-1] = 0

vorFT = np.fft.fft2(vo_en)
#vor   = np.fft.ifft2(inv_v)

kx = np.zeros(vo_en.shape[0])
ky = np.zeros(vo_en.shape[1]) 

kx = np.fft.fftfreq(vo_en.shape[0],dx/(2*np.pi))
ky = np.fft.fftfreq(vo_en.shape[1],dx/(2*np.pi))

M = np.add.outer(kx**2,ky**2)

psiFT = -1/M * vorFT

psi = np.fft.ifft2(psiFT)

u_comp = -ddy(psi,dy)
v_comp = ddx(psi,dx)


### correlation ####

#lat_b_i_u = np.where(lats < 60)[0][-1]
#lat_b_i_l = np.where(lats > 30)[0][0]
#a = vorticity2[:,lat_b_i_l:lat_b_i_u+1]
#b = vo_an[:,lat_b_i_l:lat_b_i_u+1]
#c = vorticity[:,lat_b_i_l:lat_b_i_u+1]

#fc_corr,p = stats.pearsonr(np.reshape(a,(a.shape[0]*a.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
#print("fc_corr=" + str(fc_corr))

#pers_corr, pp = stats.pearsonr(np.reshape(c,(c.shape[0]*c.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
#print("pers_corr=" + str(pers_corr))    
    
#########plotting###########    

map = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map.contourf(coord[0],coord[1],u.transpose(),np.arange(-10,10,0.2))
    
plt.tight_layout()
    
plt.savefig("u_an.png")

map2 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map2.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map2.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map2.contourf(coord[0],coord[1],v.transpose(),np.arange(-10,10,0.2))
    
plt.tight_layout()
    
plt.savefig("v_an.png")

map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map3.contourf(coord[0],coord[1],u_comp.transpose(),np.arange(-10,10,0.2))
    
plt.tight_layout()
    
plt.savefig("u_comp.png")

map4 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map4.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map4.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map4.contourf(coord[0],coord[1],v_comp.transpose(),np.arange(-10,10,0.2))
    
plt.tight_layout()
    
plt.savefig("v_comp.png")

