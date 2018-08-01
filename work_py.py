import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import matplotlib.colors as mcolors 

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
        'size'   : 14,  
       }  
######### Parameter you should set ###########
start   =  849  # start time
stop    =  849  # end time
step    =  1  # the interval or step

directory = './tacc_worky_txt/'
px_y = np.loadtxt(directory+'px2d.txt')
py_y = np.loadtxt(directory+'py2d.txt')
xx_y = np.loadtxt(directory+'xx2d.txt')
yy_y = np.loadtxt(directory+'yy2d.txt')
workx2d_y = np.loadtxt(directory+'workx2d.txt')
worky2d_y = np.loadtxt(directory+'worky2d.txt')
fieldex_y = np.loadtxt(directory+'field2dex.txt')
fieldey_y = np.loadtxt(directory+'field2dey.txt')
fieldbz_y = np.loadtxt(directory+'field2dbz.txt')
gg_y = (px_y**2+py_y**2+1)**0.5


directory = './tacc_workx_txt/'
px_x = np.loadtxt(directory+'px2d.txt')
py_x = np.loadtxt(directory+'py2d.txt')
xx_x = np.loadtxt(directory+'xx2d.txt')
yy_x = np.loadtxt(directory+'yy2d.txt')
workx2d_x = np.loadtxt(directory+'workx2d.txt')
worky2d_x = np.loadtxt(directory+'worky2d.txt')
fieldex_x = np.loadtxt(directory+'field2dex.txt')
fieldey_x = np.loadtxt(directory+'field2dey.txt')
fieldbz_x = np.loadtxt(directory+'field2dbz.txt')
gg_x = (px_x**2+py_x**2+1)**0.5


color_y = np.zeros_like(px_y[:,0])
color_x = np.zeros_like(px_x[:,0])

for i in range(0,color_y.size):
    color_y[i] = py_y[i,np.argmax(xx_y[0,:] > 20)]
for i in range(0,color_x.size):
    color_x[i] = py_x[i,np.argmax(xx_x[0,:] > 20)]



for n in range(start,stop+step,step):
    #### header data ####
#    plt.subplot()
    plt.scatter(workx2d_y[:,n], worky2d_y[:,n], c=abs(color_y), norm=colors.Normalize(vmin=0, vmax=25), s=50, cmap='rainbow', edgecolors='black', alpha=0.6)
    plt.scatter(workx2d_x[:,n], worky2d_x[:,n], c=abs(color_x), norm=colors.Normalize(vmin=0, vmax=25), s=50, cmap='rainbow', edgecolors='black', alpha=0.6)
    cbar=plt.colorbar()
    cbar.set_label(r'$P_y$ for injecting time',fontdict=font)

    plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=1.5)
    plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=1.5)
    plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),':k',linewidth=1.5)
 #   plt.legend(loc='upper right')
    plt.xlim(-200,600)
    plt.ylim(-200,600)
    plt.xlabel(r'$Work_x [m_ec^2]$',fontdict=font)
    plt.ylabel(r'$Work_y [m_ec^2]$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    #plt.show()
    #lt.figure(figsize=(100,100))
    fig = plt.gcf()
    fig.set_size_inches(8, 6.5)
    fig.savefig('./work_py/'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")

    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
