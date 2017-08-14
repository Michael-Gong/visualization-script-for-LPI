import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print 'electric field unit: '+str(exunit)
print 'magnetic field unit: '+str(bxunit)
print 'density unit nc: '+str(denunit)

font = {'family' : 'Helvetic',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

for i in [16]:
    data = sdf.read("./Datan"+str(i)+"dp/0029.sdf",dict=True)
    header=data['Header']
    time=header['time']
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    plt.subplot(1,2,1)
    name='electron'
    theta_en_Z1 = data['dist_fn/theta_en/'+name].data[:,:,0]
    theta_en_Z1=np.log(theta_en_Z1+1.0)
    #levels = np.logspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 40, base=10)
    levels = np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 40)
    dist_x  = data['Grid/theta_en/'+name].data[0]
    dist_y  = data['Grid/theta_en/'+name].data[1]/(q0*1.0e6)
    dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
    plt.contourf(dist_X, dist_Y, theta_en_Z1.T, levels=levels, cmap=cm.gnuplot2)
    #### manifesting colorbar, changing label and axis properties ####
    plt.xlim(-3.14,3.14)
    plt.ylim(0.0,120)
    plt.text(-2.55,108,r'(a) Electron',fontsize=20,color='w')
    
    
    cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
    cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
    plt.xlabel(r"$\displaystyle\theta$ [rad]",fontdict=font)
    plt.ylabel('E [MeV]',fontdict=font)
    plt.xticks(fontsize=16); plt.yticks(fontsize=16);
    plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #plt1 = plt.twinx()
    #plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)
    
    plt.subplot(1,2,2)
    name='photon'
    theta_en_Z1 = data['dist_fn/theta_en/'+name].data[:,:,0]
    theta_en_Z1=np.log(theta_en_Z1+1.0)
    #levels = np.logspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 40, base=10)
    levels = np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 40)
    dist_x  = data['Grid/theta_en/'+name].data[0]
    dist_y  = data['Grid/theta_en/'+name].data[1]/(q0*1.0e6)
    dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
    plt.contourf(dist_X, dist_Y, theta_en_Z1.T, levels=levels, cmap=cm.gnuplot2)
    #### manifesting colorbar, changing label and axis properties ####
    plt.xlim(-3.14,3.14)
    plt.ylim(0.1,15)
    plt.text(-2.55,13.5,r'(b) Photon',fontsize=20,color='w')
    cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
    cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
    plt.xlabel(r"$\displaystyle\theta$ [rad]",fontdict=font)
    plt.ylabel('E [MeV]',fontdict=font)
    plt.xticks(fontsize=16); plt.yticks(fontsize=16);
    plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    #plt1 = plt.twinx()
    #plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)
    
    fig = plt.gcf()
    fig.set_size_inches(15, 6)
    fig.savefig('./test'+str(i)+'.png',format='png',dpi=100)
    plt.close("all")
