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

n=47

data1 = sdf.read("./0047.sdf",dict=True)
header=data1['Header']
time=header['time']


#plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)

plt.subplot(2,1,1)
name='electron'
en_Z1 = data1['dist_fn/en/'+name].data[:,0,0]
dist_x1  = data1['Grid/en/'+name].data[0]/(q0*1.0e6)

plt.plot(dist_x1, en_Z1, color='blue', linewidth=2)
#plt.legend(loc='upper right',framealpha=1.0,markerscale=4.0,fontsize=20.0)
#### manifesting colorbar, changing label and axis properties ####
#plt.xlim(0.0,125)
#plt.ylim(pow(10,8),pow(10,18))
#plt.text(8,pow(10,17),r'(a) Electron',fontsize=20,color='k')
plt.yscale('log')
plt.xlabel("Energy [MeV]",fontdict=font)
plt.ylabel("dN/dE [A.U.]",fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)

plt.subplot(2,1,2)
name='photon'
en_Z1 = data1['dist_fn/en/'+name].data[:,0,0]
dist_x1  = data1['Grid/en/'+name].data[0]/(q0*1.0e6)

plt.plot(dist_x1, en_Z1, color='red', linewidth=2)
#### manifesting colorbar, changing label and axis properties ####
plt.xlim(0.0,11.92)
plt.yscale('log')
plt.xlabel("Energy [MeV]",fontdict=font)
plt.ylabel("dN/dE [A.U.]",fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)

fig = plt.gcf()
fig.set_size_inches(9, 15)
fig.savefig('./en_47.png',format='png',dpi=160)
plt.close("all")
