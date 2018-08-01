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

tt = np.linspace(5.0,89.9,850)
#tt = tt[:,np.newaxis]

select_x = np.array([21,17,2,57,58,68,72,78,172,107])


    
line_y1 = np.linspace(-12,-3.2,1001)
line_x1 = np.zeros_like(line_y1)+5
    
line_y2 = np.linspace(3.2,12,1001)
line_x2 = np.zeros_like(line_y2)+5
    
line_x3 = np.linspace(5,200.0,1001)
line_y3 = np.zeros_like(line_x3)+3.2

line_x4 = np.linspace(5,200.0,1001)
line_y4 = np.zeros_like(line_x4)-3.2

for index in select_x:
    norm_x = matplotlib.colors.Normalize(vmin=1.,vmax=501.)
    norm_y = matplotlib.colors.Normalize(vmin=1.,vmax=301.)

    plt.scatter(xx[index,50:], yy[index,50:], c=gg[index,50:], norm=norm_x, s=6, cmap='autumn', edgecolors='None')
    plt.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')


cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(501), 5) )
cbar.set_label('$\gamma$',fontdict=font)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
#   plt.legend(loc='upper right')




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

number=400

tt = np.linspace(5.0,89.9,850)
#tt = tt[:,np.newaxis]

select_x = np.array([0,1,2,3,12,25,26,29,36,28])


line_y1 = np.linspace(-12,-3.2,1001)
line_x1 = np.zeros_like(line_y1)+5
    
line_y2 = np.linspace(3.2,12,1001)
line_x2 = np.zeros_like(line_y2)+5
    
line_x3 = np.linspace(5,200.0,1001)
line_y3 = np.zeros_like(line_x3)+3.2

line_x4 = np.linspace(5,200.0,1001)
line_y4 = np.zeros_like(line_x4)-3.2

for index in select_x:
    norm_x = matplotlib.colors.Normalize(vmin=1.,vmax=501.)
    norm_y = matplotlib.colors.Normalize(vmin=1.,vmax=301.)

    plt.scatter(xx[index,50:], yy[index,50:], c=gg[index,50:], norm=norm_y, s=6, cmap='winter', edgecolors='None')
    plt.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
    plt.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')


cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(301), 5) )
cbar.set_label('$\gamma$',fontdict=font)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
#   plt.legend(loc='upper right')

plt.xlim(-5,130)
plt.ylim(-3.8,3.8)
plt.xlabel('x [$\lambda$]',fontdict=font)
plt.ylabel('y [$\lambda$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)




fig = plt.gcf()
fig.set_size_inches(16.0, 6.5)
fig.savefig('./figure_wrap_up/'+'trajectory_xy.png',format='png',dpi=160)
plt.close("all")
print('finised '+str(index))

