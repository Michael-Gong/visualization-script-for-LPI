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

data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data

#choice = np.random.choice(range(px.size), 10000, replace=False)
gamma  = work_x+work_y+1


value_axisx = np.linspace(7,700,50)
value_axisy = np.linspace(7,700,50)
value_grid = np.linspace(0,700,51)

value_total_x = np.zeros_like(value_axisy)
value_total_y = np.zeros_like(value_axisy)
value_num   = np.zeros_like(value_axisy)

for i in range(50):
    value_total_x[i] = np.sum(work_x[(value_grid[i]<=gamma) & (value_grid[i+1]>gamma)],0)
    value_total_y[i] = np.sum(work_y[(value_grid[i]<=gamma) & (value_grid[i+1]>gamma)],0)
    value_num[i] = np.size(work_y[(value_grid[i]<=gamma) & (value_grid[i+1]>gamma)])
    print('x-:',value_total_x[i]/(value_total_x[i]+value_total_y[i]),'; y-:',value_total_y[i]/(value_total_x[i]+value_total_y[i]))

#    plt.subplot()
y_x = value_total_x/(value_total_x+value_total_y)
y_x[y_x > 1] = 1
y_y = 1-y_x
width=10
pl=plt.bar(value_axisx, y_x*100, width, color='orangered',edgecolor='black',linewidth=1)
pt=plt.bar(value_axisx, y_y*100, width, bottom=y_x*100, color='dodgerblue',edgecolor='black',linewidth=1)

plt.xlim(-10,717)
plt.ylim(0,103)
plt.xlabel('$\epsilon_e$ [m$_e$c$^2$]',fontdict=font)
plt.ylabel('W$_{x(y)}$/$\epsilon_e$ [%]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.legend(['W$_x$','W$_y$'],loc='lower right',fontsize=20,framealpha=0.5)
#plt.text(200,650,' t=400fs',fontdict=font)
plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.96,
                wspace=None, hspace=None)

#plt.show()
#lt.figure(figsize=(100,100))
fig = plt.gcf()
fig.set_size_inches(7, 6.)
fig.savefig('./figure_wrap_up/work_l_t_new2'+str(n).zfill(4)+'.png',format='png',dpi=160)
#plt.close("all")

print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
