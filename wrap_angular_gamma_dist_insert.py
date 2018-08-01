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

px = px [(abs(grid_y) < 3.2) & (gg > 200.0)]
py = py [(abs(grid_y) < 3.2) & (gg > 200.0)]
work_x = work_x [(abs(grid_y) < 3.2) & (gg > 200.0)]
work_y = work_y [(abs(grid_y) < 3.2) & (gg > 200.0)]

gg = (px**2+py**2+1)**0.5
theta = np.arctan2(py,px)*180.0/np.pi

x_1 = np.linspace(-12.,12.,200)
y_1,thrush = np.histogram(theta[work_x > work_y],bins=200,range=(-12.5,12.5))
y_2,thrush = np.histogram(theta[work_x < work_y],bins=200,range=(-12.5,12.5))

plt.plot(x_1,(y_1+y_2)/1000,'-k',linewidth=3,label='Total')
plt.plot(x_1,y_1/1000,'--r',linewidth=3, label='Work$_x$ > Work$_y$')
plt.plot(x_1,y_2/1000,'--b',linewidth=3, label='Work$_x$ > Work$_y$')
plt.ylabel('dN/d'+r'$\theta$'+'[A.U.]', color='k',fontdict=font)
plt.xlabel(r'$\theta$'+' [degree]', color='k', fontdict=font)
plt.legend(loc='best',fontsize=18,framealpha=0.5)


plt.yticks(fontsize=20,color='k')
plt.xticks(fontsize=20,color='k')
plt.ylim(0,0.8)
plt.xlim(-12,12)


#plt.show()
#lt.figure(figsize=(100,100))

fig = plt.gcf()
fig.set_size_inches(7, 6.4)
fig.savefig('./figure_wrap_up/'+'angular_gamma_dist_0012_insert.png',format='png',dpi=160)
plt.close("all")

