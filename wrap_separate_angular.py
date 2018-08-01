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


data = sdf.read('./Data/'+str(12).zfill(4)+".sdf",dict=True)
header=data['Header']
time=header['time']

px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
#field_ex = data['Particles/field_ex/subset_high_e/electron'].data/exunit
#field_ey = data['Particles/field_ey/subset_high_e/electron'].data/exunit
#field_bz = data['Particles/field_bz/subset_high_e/electron'].data/bxunit
gg = (px**2+py**2+1)**0.5

px = px [(abs(grid_y) < 3.2) & (gg > 0.0)]
py = py [(abs(grid_y) < 3.2) & (gg > 0.0)]
work_x = work_x [(abs(grid_y) < 3.2) & (gg > 0.0)]
work_y = work_y [(abs(grid_y) < 3.2) & (gg > 0.0)]

gg = (px**2+py**2+1)**0.5
theta = np.arctan2(py,px)*180.0/np.pi


x_1 = np.linspace(-12.,12.,200)
y_1,thrush = np.histogram(theta[work_x > work_y],bins=200,range=(-12.5,12.5))
y_2,thrush = np.histogram(theta[work_x < work_y],bins=200,range=(-12.5,12.5))


plt.subplot(3,1,1)
#    plt.subplot()
plt.grid(which='major',color='k', linestyle='--', linewidth=0.2)
plt.scatter(theta[work_x > work_y], gg[work_x > work_y], c='red',  s=1., edgecolors='None', alpha=0.5, label='W$_x$ > W$_y$')
#plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
plt.ylabel('$\gamma$',fontdict=font)
plt.legend(loc='upper right',scatterpoints=1,markerscale=8,fontsize=16,framealpha=0.)
plt.xticks(fontsize=0); plt.yticks(fontsize=20);
plt.xlim(-12,12)
plt.ylim(25,740.0)
#plt.twinx()
#plt.plot(x_1,y_1/1000,linestyle='--',color='red',linewidth=3, label='Work$_x$ > Work$_y$')
#plt.ylabel('dN/d'+r'$\theta$'+'[A.U.]', color='blue',fontdict=font)
#plt.yticks(fontsize=20,color='blue')
#plt.ylim(0,0.9)
#plt.xlim(-12,12)


plt.subplot(3,1,2)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.2)
plt.scatter(theta[work_x < work_y], gg[work_x < work_y], c='blue', s=1., edgecolors='None', alpha=0.5, label='W$_x$ < W$_y$')
plt.legend(loc='upper right',scatterpoints=1,markerscale=8,fontsize=16,framealpha=0.)
#plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
plt.ylabel('$\gamma$',fontdict=font)
plt.xticks(fontsize=0); plt.yticks(fontsize=20);
plt.xlim(-12,12)
plt.ylim(25,740.0)
#plt.twinx()
#plt.plot(x_1,y_2/1000,linestyle='--',color='blue',linewidth=3, label='Work$_x$ > Work$_y$')
#plt.ylabel('dN/d'+r'$\theta$'+'[A.U.]', color='blue',fontdict=font)
#plt.yticks(fontsize=20,color='blue')
#plt.ylim(0,0.9)
#plt.xlim(-12,12)

#    plt.ylim(0,400)
plt.subplot(3,1,3)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.2)
plt.plot(x_1,(y_1+y_2)/1000,'-k',linewidth=3,label='Total')
plt.plot(x_1,y_1/1000,linestyle='--',color='red',linewidth=3, label='W$_x$ > W$_y$')
plt.plot(x_1,y_2/1000,linestyle='--',color='blue',linewidth=3, label='W$_x$ < W$_y$')
plt.legend(loc='upper right',fontsize=16,framealpha=0.)
plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
plt.ylabel('dN/d'+r'$\theta$'+'[A.U.]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.ylim(0,0.9)
plt.xlim(-12,12)

plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,
                wspace=0.05, hspace=0.05)


fig = plt.gcf()
fig.set_size_inches(9.2, 7.2)
fig.savefig('./figure_wrap_up/'+'angular_separate.png',format='png',dpi=160)
plt.close("all")

