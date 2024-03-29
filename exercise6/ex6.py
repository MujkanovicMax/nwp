#!/usr/bin/env python
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats



def leap_lag(values,mu,mu_D,t,dt,plot=False,dx=None,dy=None):
            
    fts= np.zeros(values.shape)
    dvalues = fts*1
    
    fts[0,:] = values[0,:] - mu[0,:]*(values[1,:]-values[-1,:])
    fts[1:-1,:] = values[1:-1,:] - mu[1:-1,:]*(values[2:,:]-values[0:-2,:])
    fts[-1,:] = values[-1,:] - mu[-1,:]*(values[0,:]-values[-2,:])
    
    for j,time in enumerate(np.arange(0,t,dt)):
        dvalues[0,:] = values[0,:] - mu[0,:] *(fts[1,:]-fts[-1,:]) +2*mu_D* (values[1,:] - 2*values[0,:] + values[-1,:])
        dvalues[1:-1,:] = values[1:-1,:] - mu[1:-1,:] *(fts[2:,:]-fts[0:-2,:]) +2*mu_D* (values[2:,:] - 2*values[1:-1,:] + values[0:-2,:])
        dvalues[-1,:] = values[-1,:] - mu[-1,:] *(fts[0,:]-fts[-2,:]) +2*mu_D* (values[0,:] - 2*values[-1,:] + values[-2,:])
        values= fts*1
        fts=dvalues*1
        
        u_fc = -ddy(dvalues,dy)
        v_fc = ddx(dvalues,dx)
        vo_fc = ddx(v_fc,dx)-ddy(u_fc,dy)
        
        if plot:
            map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
            map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
            map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
            map3.contourf(coord[0],coord[1],vo_fc.transpose(),np.arange(-4e-4,4.0001e-4,1e-6))
    
            plt.tight_layout()
                
            plt.savefig("vo_fc_"+str(time)+"_min.png")
        
    return dvalues


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

def double_ddx(values,dx):
    dvalues = np.zeros(values.shape)
    dvalues[0,:]      = values[-1,:]-2*values[0,:]+values[1,:] 
    dvalues[1:-1,:]   = values[0:-2,:]-2*values[1:-1,:] + values[2:,:]
    dvalues[-1,:]     = values[-2,:] -2*values[-1,:] + values[0,:]
    
    dvalues = dvalues / dx / dx
    
    return dvalues

def double_ddy(values,dy):
    dvalues = np.zeros(values.shape)
    dvalues[:,0]      = 0*(-values[:,0]+values[:,1])/dy/dy 
    dvalues[:,1:-1]   = (values[:,0:-2] -2*values[:,1:-1] + values[:,2:])/(dy*dy)
    dvalues[:,-1]     = 0*(-values[:,-2] + values[:,-1])/dy/dy
    
    return dvalues
def inversion(vo,umean,dx,dy):

    vo_en = np.hstack((np.fliplr(vo[:,1:-1])*(-1),vo))
    vo_en[:,vo[:,1:-1].shape[1]] = 0
    vo_en[:,-1] = 0

    vorFT = np.fft.fft2(vo_en)
    #vor   = np.fft.ifft2(inv_v)


    kx = np.fft.fftfreq(vo_en.shape[0],dx/(2*np.pi))
    ky = np.fft.fftfreq(vo_en.shape[1],dy/(2*np.pi))
    

    M = np.add.outer(kx**2,ky**2)
    M[0,0] = 1

    psiFT = -1/M * vorFT

    psi = np.fft.ifft2(psiFT)
    psi = np.real(psi)


    psi = psi[:,vo[:,1:-1].shape[1]:]-umean*dy*np.arange(0,vo.shape[1])
    
    return psi



def mkfc(values,f,u,v,D,dx,dy):
    
    
    #a=u[:,30:61]
    #c=v[:,30:61]
    
    umean = np.nanmean(u)
    
    fts = values + dt * (- u * ddx(values,dx) - v * ddy(values,dy) - v * ddy(f,dy))
    psi=inversion(fts,umean,dx,dy)
    u = -ddy(psi,dy)
    v = ddx(psi,dx)
    #print(v[:,0],v[:,-1])
    
        
    #b = u[:,30:61]
    #d = v[:,30:61]

    #fc_corr,p = stats.pearsonr(np.reshape(a,(a.shape[0]*a.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
    #print("u_corr=" + str(fc_corr))

    #pers_corr, pp = stats.pearsonr(np.reshape(c,(c.shape[0]*c.shape[1],)),np.reshape(d,(d.shape[0]*d.shape[1],)))
    #print("v_corr=" + str(pers_corr))
    
    print(np.mean(u),np.mean(v), np.mean(fts), np.mean(fts**2))
    
    for j,time in enumerate(np.arange(0,t,dt)):
            
            beta = ddy(f,dy)
            
            
            dvalues = values + 2*dt * (- u * ddx(fts,dx) - v * ddy(fts,dy) - v * beta + D * (double_ddx(values,dx) + double_ddy(values,dy)))
            #print(- u * ddx(fts,dx) - v * ddy(fts,dy) - v * ddy(f,dy) )#+ D * (ddx(ddx(values,dx),dx) + ddy(ddy(values,dy),dy)))
            values= fts*1
            fts=dvalues*1
            
            psi=inversion(dvalues,umean,dx,dy)

            #print(np.nanmean(u))
            u = -ddy(psi,dy)
            v = ddx(psi,dx)
            #print(vo_fc)
            
            print(np.mean(u),np.mean(v), np.mean(fts), np.mean(fts**2))

            #map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            #llcrnrlon=-180,urcrnrlon=180,resolution='c')
            #map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
            #map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
            #map3.contourf(coord[0],coord[1],v.transpose(),100)
    
            #plt.tight_layout()
                
            #plt.savefig("v_"+str(time)+"_min.png")
           
    return dvalues


#Read reanalysis data
filename='../exercise4/eraint_2019020100.nc'
ncf = netCDF4.Dataset(filename, 'r')
filename='../exercise4/eraint_2019020200.nc'
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
vo_error =  0.1*vo_an
u_error = u*0.1

#vo += vo_error
#u  += u_error


coord = np.meshgrid(lons-180.,lats)
radius_e = 6371000
dx = radius_e*np.pi/180.* np.cos(45/180.*np.pi)
dy = radius_e*np.pi/180.
dt = 0.1*3600.
t= 8*7*6*5*4*3*2*1*15#~7tage

mu=u*dt/dx

D = 5*10000#/np.cos(np.pi/4)/6374000/2/np.pi*360/np.cos(np.pi/4)/6374000/2/np.pi*360
mu_D=dt/dx/dx * D
omega=2*np.pi/24./3600

f = 2*omega*np.sin(lats[ind_x]/180.*np.pi)
f = np.outer(np.ones(lons.shape),f)

vo = ddx(v,dx) - ddy(u,dy)

map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map3.contourf(coord[0],coord[1],vo.transpose(),np.linspace(-4e-4,4e-4,1000))

plt.tight_layout()
    
plt.savefig("votest1.png")
plt.close('all')





#psi_fc = leap_lag(psi,mu,mu_D,t,dt,plot=True,dx=dx,dy=dy)
vo_fc = mkfc(vo,f,u,v,D,dx,dy)

#u_fc = -ddy(psi_fc,dy)
#v_fc = ddx(psi_fc,dx)

map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map3.contourf(coord[0],coord[1],vo_fc.transpose(),np.linspace(-4e-4,4e-4,1000))

plt.tight_layout()
    
plt.savefig("vofc.png")
plt.close('all')


### correlation ####

#ind = np.where((lats>=30)&(lats <= 60))
#a = u[:,ind[0]]
#b = u_comp[:,ind[0]]
##c = vorticity[:,lat_b_i_l:lat_b_i_u+1]

#fc_corr,p = stats.pearsonr(np.reshape(a,(a.shape[0]*a.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
#print("fc_corr=" + str(fc_corr))

#pers_corr, pp = stats.pearsonr(np.reshape(c,(c.shape[0]*c.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
#print("pers_corr=" + str(pers_corr))    
    
########plotting###########    
    
#map = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            #llcrnrlon=-180,urcrnrlon=180,resolution='c')
#map.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
#map.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
#map.contourf(coord[0],coord[1],u.transpose(),np.arange(-100,100,0.2))
    
#plt.tight_layout()
    
#plt.savefig("u_an.png")

map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map3.contourf(coord[0],coord[1],vo_fc.transpose(),np.linspace(-4e-4,4e-4,1000))#,np.arange(-4e-6,4.0001e-4,1e-6))

plt.tight_layout()
    
plt.savefig("vo_fc_"+str(time)+"_min.png")
plt.close('all')

map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
map3.contourf(coord[0],coord[1],vo_an.transpose(),np.linspace(-4e-4,4e-4,1000))#,np.arange(-4e-6,4.0001e-4,1e-6))

plt.tight_layout()
    
plt.savefig("vo_an"+".png")
plt.close('all')

#map3 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            #llcrnrlon=-180,urcrnrlon=180,resolution='c')
#map3.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
#map3.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
#map3.contourf(coord[0],coord[1],u_fc.transpose(),np.arange(-100,100,0.2))
    
#plt.tight_layout()
    
#plt.savefig("u_fc.png")

#map2 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            #llcrnrlon=-180,urcrnrlon=180,resolution='c')
#map2.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
#map2.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
#map2.contourf(coord[0],coord[1],v.transpose(),np.arange(-100,100,0.2))
    
#plt.tight_layout()
    
#plt.savefig("v_an.png")

#map4 = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            #llcrnrlon=-180,urcrnrlon=180,resolution='c')
#map4.drawmapboundary(color='k', linewidth=1.0, fill_color=None)
#map4.drawcoastlines(linewidth=1.0, linestyle='solid', color='k')
#map4.contourf(coord[0],coord[1],v_fc.transpose(),np.arange(-100,100,0.2))
    
#plt.tight_layout()
    
#plt.savefig("v_fc.png")

a=vo_fc[:,30:61]
b=vo_an[:,30:61]

fc_corr,p = stats.pearsonr(np.reshape(a,(a.shape[0]*a.shape[1],)),np.reshape(b,(b.shape[0]*b.shape[1],)))
print("u_corr=" + str(fc_corr))

#pers_corr, pp = stats.pearsonr(np.reshape(c,(c.shape[0]*c.shape[1],)),np.reshape(d,(d.shape[0]*d.shape[1],)))
#print("v_corr=" + str(pers_corr))

#print(np.mean(u),np.mean(v), np.mean(fts), np.mean(fts**2))






