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
        'size'   : 30,  
       }  
######### Parameter you should set ###########
start   =  1  # start time
stop    =  49  # end time
step    =  1  # the interval or step

n=12

px = np.loadtxt('./txt/px_'+str(n).zfill(4)+'sdf.txt')
py = np.loadtxt('./txt/py_'+str(n).zfill(4)+'sdf.txt')
grid_x = np.loadtxt('./txt/grid_x_'+str(n).zfill(4)+'sdf.txt')
grid_y = np.loadtxt('./txt/grid_y_'+str(n).zfill(4)+'sdf.txt')
work_x = np.loadtxt('./txt/work_x_'+str(n).zfill(4)+'sdf.txt')
work_y = np.loadtxt('./txt/work_y_'+str(n).zfill(4)+'sdf.txt')

#choice = np.random.choice(range(px.size), 10000, replace=False)
choice = np.random.choice(range(px.size), px.size, replace=False)
px = px[choice]
py = py[choice]
work_x = work_x[choice]
work_y = work_y[choice]


theta = np.arctan2(py,px)*180.0/np.pi
    

theta[theta < -7.5] = -7.5
theta[theta >  7.5] =  7.5


color_index = abs(theta)

#    plt.subplot()
plt.scatter(work_x, work_y, c=color_index, s=3, cmap='rainbow_r', edgecolors='None', alpha=0.66)
cbar=plt.colorbar( ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.01)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
cbar.set_label(r'$|\theta|$'+' [degree]',fontdict=font)

plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
 #   plt.legend(loc='upper right')
plt.xlim(-250,750)
plt.ylim(-250,750)
plt.xlabel('W$_x$ [m$_e$c$^2$]',fontdict=font)
plt.ylabel('W$_y$ [m$_e$c$^2$]',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
plt.text(-100,650,' t = 400 fs',fontdict=font)
plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                wspace=None, hspace=None)

#plt.show()
#lt.figure(figsize=(100,100))
fig = plt.gcf()
fig.set_size_inches(12, 10.5)
fig.savefig('./figure_wrap_up/theta_new'+str(n).zfill(4)+'.png',format='png',dpi=160)
#plt.close("all")

print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
