#%matplotlib inline
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
#from colour import Color

######## Constant defined here ########
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

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  



part_number=10000
nsteps=2
axis_a=np.linspace(10,160,16)
axis_w=np.linspace(30,90,16)
rebo1=np.zeros([16,16])
rebo2=np.zeros([16,16])

insert1='epoch2drr/'
insert2='epoch2dqe/'
n=12
for ia in range(16):
    for iw in range(16):
	data  = sdf.read(insert1+'Dataa'+str(int(axis_a[ia]))+'w'+str(int(axis_w[iw]))+'/'+str(n).zfill(4)+'.sdf',dict=True)
    	px    = data['Particles/Px/electron'].data/(m0*v0)
    	py    = data['Particles/Py/electron'].data/(m0*v0)
        rebo1[ia,iw]=float(np.size(py[py > 0.50]))/np.size(py)
	data  = sdf.read(insert2+'Dataa'+str(int(axis_a[ia]))+'w'+str(int(axis_w[iw]))+'/'+str(n).zfill(4)+'.sdf',dict=True)
    	px    = data['Particles/Px/electron'].data/(m0*v0)
    	py    = data['Particles/Py/electron'].data/(m0*v0)
        rebo2[ia,iw]=float(np.size(py[py > 0.50]))/np.size(py)

plt.subplot(1,2,1) 
X,Y=np.meshgrid(axis_a,axis_w/10.)
Z=rebo1
#plt.imshow(rebo1,interpolation='none', cmap=cm.nipy_spectral)
levels = np.linspace(0, 0.5, 40)
plt.contourf(X, Y, Z.T, levels=levels, cmap=cm.terrain)
#### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=[0,0.25,0.5])
cbar.set_label('Reflection ratio',fontdict=font)
plt.xlabel('$a_0$',fontdict=font)
plt.ylabel('$r_0 [\mu m]$',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('LL RR scan')


plt.subplot(1,2,2) 
X,Y=np.meshgrid(axis_a,axis_w/10.)
Z=rebo2
#plt.imshow(rebo1,interpolation='none', cmap=cm.nipy_spectral)
levels = np.linspace(0, 0.5, 40)
plt.contourf(X, Y, Z.T, levels=levels, cmap=cm.terrain)
#### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=[0,0.25,0.5])
cbar.set_label('Reflection ratio',fontdict=font)
plt.xlabel('$a_0$',fontdict=font)
plt.ylabel('$r_0 [\mu m]$',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('QED RR scan')


fig = plt.gcf()
fig.set_size_inches(20, 8)
fig.savefig('scan_a_w.png',format='png',dpi=100)
plt.close("all")


