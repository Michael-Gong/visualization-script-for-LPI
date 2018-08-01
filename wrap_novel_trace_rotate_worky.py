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

directory = './txt_1300/'
px = np.loadtxt(directory+'px2d_x.txt')
py = np.loadtxt(directory+'py2d_x.txt')
xx = np.loadtxt(directory+'xx2d_x.txt')
yy = np.loadtxt(directory+'yy2d_x.txt')
workx2d = np.loadtxt(directory+'workx2d_x.txt')
worky2d = np.loadtxt(directory+'worky2d_x.txt')
fieldex = np.loadtxt(directory+'fieldex2d_x.txt')/4.0
fieldey = np.loadtxt(directory+'fieldey2d_x.txt')/4.0
fieldbz = np.loadtxt(directory+'fieldbz2d_x.txt')/4.0

ey_averaged = -8.0/3.2*yy
bz_averaged = -8.0/3.2*yy

laser_ey = fieldey-ey_averaged
laser_bz = fieldbz-bz_averaged


gg = (px**2+py**2+1)**0.5
R = gg-px
theta = np.arctan2(py,px)

number=400



plt.subplot(2,1,1) 
#axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left')
#axin2 = inset_axes(ax, width='15%', height='5%', loc='upper center')

grid_t = np.zeros(130)
grid_y = np.linspace(-3.5,3.5,140)
grid_data = np.zeros([130,140])
index=68
for n in range(5,130,1):
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    xi=np.min(np.where(xx[index,(n-5)*10] < x))
    print(np.shape(ex[xi, 380:520]))
    grid_data[n,:] = (ex[xi, 380:520] + ex[xi-1, 380:520])*0.5


X, Y = np.meshgrid(grid_t, grid_y) 
levels = np.linspace(-28.1, 28.1, 50)
grid_data[grid_data < -28]=-28
grid_data[grid_data >  28]= 28
plt.contourf((Y).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-28,28,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)


norm_x = matplotlib.colors.Normalize(vmin=np.min(worky2d[index,0:1299:20]),vmax=np.max(worky2d[index,0:1299:20]))
print(np.shape(grid_t),np.shape(xx[index,0:1299:10]-(n-15)))
plt.scatter(yy[index,0:1299:20], grid_t[5::2], c=worky2d[index,0:1299:20], norm=norm_x, s=60, cmap='hot', edgecolors='black')
cbar=plt.colorbar()#orientation="horizontal")
cbar.set_label('Work$_y$ [$m_ec^2$]', fontdict=font2)

#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('y [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20)
plt.text(1,350,'LDA',fontdict=font)
#plt.yscale('log')
plt.xlim(-3.4,3.4)
plt.ylim(10,430)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)










directory = './txt_1300/'
px = np.loadtxt(directory+'px2d_y.txt')
py = np.loadtxt(directory+'py2d_y.txt')
xx = np.loadtxt(directory+'xx2d_y.txt')
yy = np.loadtxt(directory+'yy2d_y.txt')
workx2d = np.loadtxt(directory+'workx2d_y.txt')
worky2d = np.loadtxt(directory+'worky2d_y.txt')
fieldex = np.loadtxt(directory+'fieldex2d_y.txt')/4.0
fieldey = np.loadtxt(directory+'fieldey2d_y.txt')/4.0
fieldbz = np.loadtxt(directory+'fieldbz2d_y.txt')/4.0

ey_averaged = -8.0/3.2*yy
bz_averaged = -8.0/3.2*yy

laser_ey = fieldey-ey_averaged
laser_bz = fieldbz-bz_averaged


gg = (px**2+py**2+1)**0.5
R = gg-px
theta = np.arctan2(py,px)


plt.subplot(2,1,2) 

grid_t = np.zeros(130)
grid_y = np.linspace(-3.5,3.5,140)
grid_data = np.zeros([130,140])
index=3
for n in range(5,130,1): 
    data = sdf.read("./Data_a20_130_2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    grid_t[n]=header['time']/1e-15

    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    xi=np.min(np.where(xx[index,(n-5)*10] < x))
    print(np.shape(ex[xi, 380:520]))
    grid_data[n,:] = (ex[xi, 380:520] + ex[xi-1, 380:520])*0.5


X, Y = np.meshgrid(grid_t, grid_y) 
levels = np.linspace(-28.1, 28.1, 50)
grid_data[grid_data < -28]=-28
grid_data[grid_data >  28]= 28
plt.contourf((Y).T, X.T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.jet)
cbar=plt.colorbar(ticks=np.linspace(-28,28,5))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega_0/e$]', fontdict=font2)
norm_x = matplotlib.colors.Normalize(vmin=np.min(worky2d[index,0:1299:20]),vmax=np.max(worky2d[index,0:1299:20]))
plt.scatter(yy[index,0:1299:20], grid_t[5::2], c=worky2d[index,0:1299:20], norm=norm_x, s=60, cmap='hot', edgecolors='black')
cbar=plt.colorbar()#orientation="horizontal")
cbar.set_label('Work$_y$ [$m_ec^2$]', fontdict=font2)
#plt.scatter(grid_t[5:], xx[index,0:1299:10]-(grid_t[5:]/3.3333333-15), c=workx2d[index,0:1299:10], norm=norm_x, s=30, cmap='hot', edgecolors='black')
#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('y [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.xlim(-3.4,3.4)
plt.ylim(10,430)
plt.text(1,350,'TDA',fontdict=font)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)


fig = plt.gcf()
fig.set_size_inches(11.2, 13.3)
fig.savefig('./figure_wrap_up/'+'move_window_rotate_ey_worky.png',format='png',dpi=160)
plt.close("all")

