import numpy as np
import matplotlib.pyplot as plt

def euler_up(phi_num, mu,t,dt):
    
    for time in np.arange(0,t,dt):
        for pos in reversed(range(len(phi_num))):
            phi_num[pos] = phi_num[pos] - mu[pos] * (phi_num[pos]-phi_num[pos-1])
    
def euler_dn(phi_num, mu,t,dt):
    
    for time in np.arange(0,t,dt):
        for pos in range(len(phi_num)):
            if pos >= 359:
                pos = pos - 360
            phi_num[pos] = phi_num[pos] - mu * (-phi_num[pos]+phi_num[pos+1])


def euler_cen(phi_num,mu,t,dt):
    for time in np.arange(0,t,dt):
        for pos in range(len(phi_num)):
            if pos >=359:
                pos=pos-360
            phi_num[pos] =phi_num[pos] - mu *(phi_num[pos+1] - phi_num[pos-1])
            
            
def leapfrog(phi_num,mu,t,dt):
    tmp=phi_num*1
    euler_up(phi_num,mu,dt,dt)
    for jj,time in enumerate(np.arange(0,t,dt)):
        for pos in range(len(phi_num)):
            if pos >=359:
                pos = pos-360
            tmp2=phi_num*1
            phi_num[pos] = tmp[pos] - mu * (tmp2[pos+1] - tmp2[pos-1])
            tmp=tmp2*1
        if(jj%50==0):
            ax.plot(x,phi_num, color = "b")
            ax.plot(x,ff,color="r")


            plt.savefig("bild_"+str(jj)+".png")
            ax.clear()

def leapfrog_diffuse_lag(phi_num,mu,mu_D,t,dt):
    tmp=phi_num*1
    euler_up(phi_num,mu,dt,dt)
    for jj,time in enumerate(np.arange(0,t,dt)):
        for pos in range(len(phi_num)):
            if pos >=359:
                pos = pos-360
            tmp2=phi_num*1
            phi_num[pos] = tmp[pos] - mu[pos] * (tmp2[pos+1] - tmp2[pos-1]) + 2*mu_D * (tmp[pos+1]-2*tmp[pos]+tmp[pos-1])
            tmp=tmp2*1
        if(jj%1==0):
            ax.plot(x,phi_num, color = "b")
            ax.plot(x,ff,color="r")


            plt.savefig("bild_"+str(jj)+".png")
            ax.clear()

def leapfrog_diffuse(phi_num,mu,mu_D,t,dt):
    tmp=phi_num*1
    euler_up(phi_num,mu,dt,dt)
    for jj,time in enumerate(np.arange(0,t,dt)):
        for pos in range(len(phi_num)):
            if pos >=359:
                pos = pos-360
            tmp2=phi_num*1
            phi_num[pos] = tmp[pos] - mu * (tmp2[pos+1] - tmp2[pos-1]) + 2*mu_D * (tmp2[pos+1]-2*tmp2[pos]+tmp2[pos-1])
            tmp=tmp2*1
        if(jj%5==0):
            ax.plot(x,phi_num, color = "b")
            ax.plot(x,ff,color="r")


            plt.savefig("bild_"+str(jj)+".png")
            ax.clear()



t = 10*24*3600
#u = 10/np.cos(np.pi/4)/6374000/2/np.pi*360


sigma = 10
phi0 = 500
x0 = 60
x = np.linspace(0,359,num=360)
u = 10 *np.exp(-(x-x0)*(x-x0)/(2*5*5))/np.cos(np.pi/4)/6374000/2/np.pi*360
phi_0=phi0 * np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
phi = phi0 * np.exp(-((x-u*t)-x0)*((x-u*t)-x0)/(2*sigma*sigma))





phi_num = 0.01*np.random.rand(360)*phi0 * np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))

#phi_num = phi0* np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
ff = u*1#phi_num*1

dx = 1
dt = dx/np.nanmax(u)*0.2
mu = u*dt/dx
D = 0#100000/np.cos(np.pi/4)/6374000/2/np.pi*360/np.cos(np.pi/4)/6374000/2/np.pi*360
mu_D=dt/dx/dx * D


  
  
fig = plt.figure()
ax = plt.axes()


leapfrog_diffuse_lag(u,mu,mu_D,t,dt)        





#plt.show()  

#print(np.max(phi_num))