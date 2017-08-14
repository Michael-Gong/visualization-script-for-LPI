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

plt.subplots_adjust(left=0.10,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)

plt.subplot(2,2,1)
data = sdf.read("./Datan16dp/0018.sdf",dict=True)
header=data['Header']
time=header['time']
x  = data['Grid/Grid'].data[0]/1.0e-6
x  = (x[1:]+x[:-1])*0.5
y  = data['Grid/Grid'].data[1]/1.0e-6
y  = (y[1:]+y[:-1])*0.5
X,Y= np.meshgrid(x, y)
ex = data['Electric Field/Ey'].data/exunit
eee=np.max([-np.min(ex.T),np.max(ex.T)])
levels = np.linspace(-eee, eee, 40)
Z=ex.T
plt.contourf(X[:,900:2699], Y[:,900:2699], Z[:,900:2699], levels=levels, cmap=cm.RdBu)
### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
cbar.set_label(r'$Ez [m_ec\omega/e]$',fontdict=font)        
plt.xlabel(r'X [$\mu m$]',fontdict=font)
plt.ylabel(r'Y [$\mu m$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('Ey at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
plt.text(35,15,r'(a) Ey',fontsize=20,color='k')


plt.subplot(2,2,2)
data = sdf.read("./Datan16dp/0018.sdf",dict=True)
header=data['Header']
time=header['time']
x  = data['Grid/Grid'].data[0]/1.0e-6
x  = (x[1:]+x[:-1])*0.5
y  = data['Grid/Grid'].data[1]/1.0e-6
y  = (y[1:]+y[:-1])*0.5
X,Y= np.meshgrid(x, y)
den = data['Derived/Number_Density/electron'].data/denunit
levels = np.linspace(np.min(den.T), 0.3*np.max(den.T), 40) 
Z=den.T
plt.contourf(X[:,900:2699], Y[:,900:2699], Z[:,900:2699], levels=levels, cmap=cm.gnuplot2)
#### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), 0.3*np.max(den.T), 5))
cbar.set_label(r'$n_e$ [$n_c$]', fontdict=font)
plt.xlabel(r'X [$\mu m$]',fontdict=font)
plt.ylabel(r'Y [$\mu m$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('Electron at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
plt.text(35,15,r'(b) Electron density',fontsize=20,color='w')



plt.subplot(2,2,3)
data = sdf.read("./Datan16dp/0018.sdf",dict=True)
header=data['Header']
time=header['time']
x  = data['Grid/Grid'].data[0]/1.0e-6
x  = (x[1:]+x[:-1])*0.5
y  = data['Grid/Grid'].data[1]/1.0e-6
y  = (y[1:]+y[:-1])*0.5
X,Y= np.meshgrid(x, y)
ex = data['Magnetic Field/Bz_averaged'].data/bxunit
eee=np.max([-np.min(ex.T),np.max(ex.T)])
levels = np.linspace(-eee, eee, 40)
Z=ex.T
plt.contourf(X[:,900:2699], Y[:,900:2699], Z[:,900:2699], levels=levels, cmap=cm.RdBu)
### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
cbar.set_label(r'$Bz_{averaged} [m_e\omega/e]$',fontdict=font)        
plt.xlabel(r'X [$\mu m$]',fontdict=font)
plt.ylabel(r'Y [$\mu m$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('Bz averaged at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
plt.text(35,15,r'(c) $Bz_{averaged}$',fontsize=20,color='k')


plt.subplot(2,2,4)
data = sdf.read("./Datan16dp/0029.sdf",dict=True)
header=data['Header']
time=header['time']
x  = data['Grid/Grid'].data[0]/1.0e-6
x  = (x[1:]+x[:-1])*0.5
y  = data['Grid/Grid'].data[1]/1.0e-6
y  = (y[1:]+y[:-1])*0.5
X,Y= np.meshgrid(x, y)
den = data['Derived/Number_Density/photon'].data/denunit
levels = np.linspace(np.min(den.T), np.max(den.T), 40) 
plt.contourf(X, Y, den.T, levels=levels, cmap=cm.gnuplot2)
#### manifesting colorbar, changing label and axis properties ####
cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
cbar.set_label(r'$n_{photon}$ [$n_c$]', fontdict=font)
plt.xlabel(r'X [$\mu m$]',fontdict=font)
plt.ylabel(r'Y [$\mu m$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title('Photon at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
plt.text(10,15,r'(d) Photon density',fontsize=20,color='w')





fig = plt.gcf()
fig.set_size_inches(15, 10)
fig.savefig('./test.png',format='png',dpi=100)
plt.close("all")
