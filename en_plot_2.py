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

electron=np.zeros(30)
carbon=np.zeros(30)
photon=np.zeros(30)
particle=np.zeros(30)
field=np.zeros(30)

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)

for n in np.arange(30):
	data = sdf.read("./Datan2w33p/"+str(n).zfill(4)+".sdf",dict=True)
        if data.has_key('dist_fn/en/photon'):
		Z=data['dist_fn/en/photon'].data[:,0,0]
		X=data['Grid/en/photon'].data[0]
	        photon[n]=np.sum(Z*X,0)
	Z=data['dist_fn/en/electron'].data[:,0,0]
	X=data['Grid/en/electron'].data[0]
	electron[n]=np.sum(Z*X,0)
	Z=data['dist_fn/en/carbon'].data[:,0,0]
	X=data['Grid/en/carbon'].data[0]
	carbon[n]=np.sum(Z*X,0)
        field[n]=data['Total Field Energy in Simulation (J)'].data
        particle[n]=data['Total Particle Energy in Simulation (J)'].data
plt.subplot(1,2,1)
time_x=np.arange(30)*3.3333*5
plt.plot(time_x[4:], electron[4:]/field[4]*100, label=r"electron", linewidth=3)
plt.plot(time_x[4:], carbon[4:]/field[4]*100, label=r"carbon", linewidth=3)
plt.plot(time_x[4:], photon[4:]/field[4]*1000*100, label=r"photon*1000", linewidth=3)
plt.plot(time_x[4:], (field[4]-electron[4:]-carbon[4:]-photon[4:])/field[4]*100, label=r"laser", linewidth=3)
plt.legend(loc='center right',framealpha=1.0,markerscale=3.0,fontsize=14.0)
#### manifesting colorbar, changing label and axis properties ####
plt.ylim(0.0,100.0)
plt.xlim(50,500)
plt.text(100,90,r'(a) $n_e=0.2n_c$',fontsize=20,color='k')

#cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
#cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
#plt.yscale('log')
plt.xlabel("Time [fs]",fontdict=font)
plt.ylabel("Conversion efficiency [\%]",fontdict=font)
plt.xticks(fontsize=16); plt.yticks(fontsize=16);
plt.title('Energy evolution for $n_e=0.2n_c$',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)

electron=np.zeros(30)
carbon=np.zeros(30)
photon=np.zeros(30)
particle=np.zeros(30)
field=np.zeros(30)
for n in np.arange(30):
	data = sdf.read("./Datan20w33p/"+str(n).zfill(4)+".sdf",dict=True)
        if data.has_key('dist_fn/en/photon'):
		Z=data['dist_fn/en/photon'].data[:,0,0]
		X=data['Grid/en/photon'].data[0]
	        photon[n]=np.sum(Z*X,0)
	Z=data['dist_fn/en/electron'].data[:,0,0]
	X=data['Grid/en/electron'].data[0]
	electron[n]=np.sum(Z*X,0)
	Z=data['dist_fn/en/carbon'].data[:,0,0]
	X=data['Grid/en/carbon'].data[0]
	carbon[n]=np.sum(Z*X,0)
        field[n]=data['Total Field Energy in Simulation (J)'].data
        particle[n]=data['Total Particle Energy in Simulation (J)'].data
plt.subplot(1,2,2)
time_x=np.arange(30)*3.3333*5
plt.plot(time_x[4:], (field[4]-field[4:]-carbon[4:])/field[4]*100.0, label=r"electron", linewidth=3)
plt.plot(time_x[4:], carbon[4:]/field[4]*100.0, label=r"carbon", linewidth=3)
plt.plot(time_x[4:], photon[4:]/field[4]*1000*100.0, label=r"photon*1000", linewidth=3)
#plt.plot(time_x[3:], (field[3]-electron[3:]-carbon[3:]-photon[3:])/field[3], label=r"laser", linewidth=3)
plt.plot(time_x[4:], (field[4:])/field[4]*100.0, label=r"laser", linewidth=3)
plt.legend(loc='center right',framealpha=1.0,markerscale=3.0,fontsize=14.0)
#### manifesting colorbar, changing label and axis properties ####
plt.ylim(0.0,100)
plt.xlim(50,500)
plt.text(100,90,r'(b) $n_e=2.0n_c$',fontsize=20,color='k')

#cbar=plt.colorbar(ticks=np.linspace(np.min(theta_en_Z1.T), np.max(theta_en_Z1.T), 5))
#cbar.set_label(r"$\displaystyle log_{10}dN/d\theta dE$ [A.U.]", fontdict=font)
#plt.yscale('log')
plt.xlabel("Time [fs]",fontdict=font)
plt.ylabel("Conversion efficiency [\%]",fontdict=font)
plt.xticks(fontsize=16); plt.yticks(fontsize=16);
plt.title('energy evolution for $n_e=2.0n_c$',fontdict=font)
#plt1 = plt.twinx()
#plt1.plot(dist_x,np.sum(theta_en_Z1,axis=1),'-y',linewidth=2.5)

fig = plt.gcf()
fig.set_size_inches(15, 6)
fig.savefig('./test.png',format='png',dpi=100)
plt.close("all")
