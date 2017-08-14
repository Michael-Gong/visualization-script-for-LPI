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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)

plt.subplot(1,2,1)


xxx=np.arange(7)
for n in xxx[::-1]:
	data=sdf.read("./Datan"+str(pow(2,n))+"dp/0029.sdf",dict=True)
	Z=data['dist_fn/theta_en/photon'].data[:,:,0]
	X=data['Grid/theta_en/photon'].data[0]/np.pi*180.0
        Z=np.sum(Z,1)/(X[-1]-X[-2])
	plt.plot(X,Z,label=r"$n_e="+str(pow(2,n)/100.0)+"n_c$",linewidth=2)         
plt.legend(loc='upper right',framealpha=0.50,markerscale=4.0,fontsize=14.0)
#### manifesting colorbar, changing label and axis properties ####
plt.xlim(-180,180)
plt.ylim(0.0,3.2e14)
plt.text(-180+36,0.9*3.2e14,r'(a)',fontsize=20,color='k')

#cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
#cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
#plt.yscale('log')
plt.xlabel(r"$\theta [degree]$",fontdict=font)
plt.ylabel(r"$dN/d\theta [A.U.]$",fontdict=font)
plt.xticks(fontsize=16); plt.yticks(fontsize=16);
#plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)

data=sdf.read("./Datan16dp/0003.sdf",dict=True)
total=data['Total Field Energy in Simulation (J)'].data

density=np.array([0.01,0.02,0.04,0.08,0.16,0.32,0.64])
photon1=np.zeros(7)
photon2=np.zeros(7)
photon3=np.zeros(7)
for n in np.arange(7):
	data=sdf.read("./Datan"+str(pow(2,n))+"dp/0029.sdf",dict=True)
	Z=data['dist_fn/en/photon'].data[:,0,0]
	X=data['Grid/en/photon'].data[0]
	photon1[n]=np.sum(Z*X,0)/total*100
	photon2[n]=np.sum(Z[X > 1.0*q0*1e6]*X[X > 1.0*q0*1e6],0)/total*100
	photon3[n]=np.sum(Z[X > 10*q0*1e6]*X[X > 10*q0*1e6],0)/total*100*20

plt.subplot(1,2,2)
plt.plot(density, photon1, 'o-', markersize=12, label=r"$\gamma_{ph}>0.1MeV$", linewidth=2)
plt.plot(density, photon2, 'o-', markersize=12, label=r"$\gamma_{ph}>1.0MeV$", linewidth=2)
plt.plot(density, photon3, 'o-', markersize=12, label=r"$\gamma_{ph}>10.0MeV *20$", linewidth=2)
plt.legend(loc='upper right',framealpha=0.50,markerscale=0.5,fontsize=14.0)
#### manifesting colorbar, changing label and axis properties ####
plt.xlim(0.0,0.66)
plt.ylim(0.00,0.139)
plt.text(0.66*0.1,0.139*0.9,r'(b)',fontsize=20,color='k')

#cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
#cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
#plt.xscale('log')
plt.xlabel(r"$n_e [n_c]$",fontdict=font)
plt.ylabel(r"Coversion efficiency [\%]",fontdict=font)
plt.xticks(fontsize=16); plt.yticks(fontsize=16);
#plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)


fig = plt.gcf()
fig.set_size_inches(13.3, 6)
fig.savefig('./test.png',format='png',dpi=100)
plt.close("all")
