import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage



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
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  

font2 = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }


plt.subplot(2,2,1) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
for n in range(5,130,1):
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ey = data['Electric Field/'+str.capitalize(name)].data/exunit
    ey_ave = data['Electric Field/'+str.capitalize('ey_averaged')].data/exunit
    ex = ey-ey_ave
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 449] + ex[0+(n-5)*30:150+(n-5)*30,450])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-22.1, 22.1, 50)
grid_data[grid_data < -22]=-22
grid_data[grid_data >  22]= 22
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-20,20,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=0$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,2) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
for n in range(5,130,1):
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ey = data['Electric Field/'+str.capitalize(name)].data/exunit
    ey_ave = data['Electric Field/'+str.capitalize('ey_averaged')].data/exunit
    ex = ey-ey_ave
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 469] + ex[0+(n-5)*30:150+(n-5)*30,470])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-17.1, 17.1, 50)
grid_data[grid_data < -17]=-17
grid_data[grid_data >  17]= 17
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-16,16,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=1$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,3) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
index=3
for n in range(5,130,1): 
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ey = data['Electric Field/'+str.capitalize(name)].data/exunit
    ey_ave = data['Electric Field/'+str.capitalize('ey_averaged')].data/exunit
    ex = ey-ey_ave
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 489] + ex[0+(n-5)*30:150+(n-5)*30,490])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-12.1, 12.1, 50)
grid_data[grid_data < -12]=-12
grid_data[grid_data >  12]= 12
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-12,12,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=2$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,4) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
index=3
for n in range(5,130,1): 
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ey = data['Electric Field/'+str.capitalize(name)].data/exunit
    ey_ave = data['Electric Field/'+str.capitalize('ey_averaged')].data/exunit
    ex = ey-ey_ave
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 489] + ex[0+(n-5)*30:150+(n-5)*30,490])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-6.1, 6.1, 50)
grid_data[grid_data < -6]=-6
grid_data[grid_data >  6]= 6
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-6,6,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=3$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

fig = plt.gcf()
fig.set_size_inches(24, 27)
fig.savefig('./jpg_v_ph/'+'move_window_ey.png',format='png',dpi=160)
plt.close("all")





plt.subplot(2,2,1) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
for n in range(5,130,1):
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ex'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 449] + ex[0+(n-5)*30:150+(n-5)*30,450])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-2.1, 2.1, 50)
grid_data[grid_data < -2]=-2
grid_data[grid_data >  2]= 2
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-2,2,5))#,orientation="horizontal")
cbar.set_label('$E_x$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=0$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,2) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
for n in range(5,130,1):
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ex'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 469] + ex[0+(n-5)*30:150+(n-5)*30,470])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-2.1, 2.1, 50)
grid_data[grid_data < -2]=-2
grid_data[grid_data >  2]= 2
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-2,2,5))#,orientation="horizontal")
cbar.set_label('$E_x$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=1$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,3) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
index=3
for n in range(5,130,1): 
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ex'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 489] + ex[0+(n-5)*30:150+(n-5)*30,490])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-2.1, 2.1, 50)
grid_data[grid_data < -2]=-2
grid_data[grid_data >  2]= 2
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-2,2,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=2$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

plt.subplot(2,2,4) 
grid_t = np.zeros(130)
grid_x = np.linspace(10,14.9833,150)
grid_data = np.zeros([130,150])
index=3
for n in range(5,130,1): 
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ey = data['Electric Field/'+str.capitalize(name)].data/exunit
    ey_ave = data['Electric Field/'+str.capitalize('ey_averaged')].data/exunit
    ex = ey-ey_ave
    grid_data[n,:] = (ex[0+(n-5)*30:150+(n-5)*30, 489] + ex[0+(n-5)*30:150+(n-5)*30,490])*0.5
X, Y = np.meshgrid(grid_t, grid_x) 
levels = np.linspace(-2.1, 2.1, 50)
grid_data[grid_data < -2]=-2
grid_data[grid_data >  2]= 2
plt.contourf((Y-10).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-2,2,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta_x$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(0,5)
plt.ylim(60,430)
plt.title('Cross at y=3$\mu m$',fontdict=font)
#plt.text(1,350,'LDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)

fig = plt.gcf()
fig.set_size_inches(24, 27)
fig.savefig('./jpg_v_ph/'+'move_window_ex.png',format='png',dpi=160)
plt.close("all")

