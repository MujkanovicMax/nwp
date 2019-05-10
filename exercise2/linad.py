import numpy as np
import matplotlib.pyplot as plt

def euler_up(phi_num, mu,t,dt):
    
    for time in np.arange(0,t,dt):
        for pos in reversed(range(len(phi_num))):
            phi_num[pos] = phi_num[pos] - mu * (phi_num[pos]-phi_num[pos-1])
    
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
    for time in np.arange(0,t,dt):
        for pos in range(len(phi_num)):
            if pos >=359:
                pos = pos-360
            tmp2=phi_num*1
            phi_num[pos] = tmp[pos] - mu * (tmp2[pos+1] - tmp2[pos-1])
            tmp=tmp2*1


#t = 10*24*3600
#u = 10/np.cos(np.pi/4)/6374000/2/np.pi*360
#sigma = 10
#phi0 = 500
#x0 = 60
#x = np.linspace(0,359,num=360)
#phi_0=phi0 * np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
#phi = phi0 * np.exp(-((x-u*t)-x0)*((x-u*t)-x0)/(2*sigma*sigma))







#phi_num = phi0 * np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))

#dx = 1
#dt = dx/u
#mu = u*dt/dx


#euler_up(phi_num,mu,t,dt)        
  
  
#fig = plt.figure()
#ax = plt.axes()

#ax.plot(x,phi_0,color="r")
#ax.plot(x,phi_num, color = "b")



#plt.show()  

#print(np.max(phi_num))