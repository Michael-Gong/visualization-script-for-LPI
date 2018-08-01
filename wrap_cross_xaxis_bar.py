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
        'size'   : 22,  
       }  
######### Parameter you should set ###########
start   =  1  # start time
stop    =  49  # end time
step    =  1  # the interval or step

n=12

directory = './txt_10july/'
xx = np.loadtxt(directory+'xx2d_x.txt')
yy = np.loadtxt(directory+'yy2d_x.txt')
workx2d = np.loadtxt(directory+'workx2d_x.txt')
worky2d = np.loadtxt(directory+'worky2d_x.txt')
gamma = workx2d+worky2d+1
#choice = np.random.choice(range(px.size), 10000, replace=False)
#choice = np.random.choice(range(xx.size), xx.size, replace=False)
#xx = xx[choice]
#yy = yy[choice]
#work_x = work_x[choice]
#work_y = work_y[choice]

cross_or_not = np.zeros_like(range(yy[:,-1].size))

for i in range(yy[:,-1].size):
    yy1 = yy[i,:-1]
    if np.size(yy1[yy[i,:-1]*yy[i,1:] < 0.0])>0:
       cross_or_not[i] = 1.0
       print(i, 'is the crossed !')
    

value_axisx = np.linspace(20,700,50)
value_axisy = np.linspace(20,700,50)
value_grid = np.linspace(20,700,51)

value_total_x_cross = np.zeros_like(value_axisy)
value_total_x_not   = np.zeros_like(value_axisy)
value_total_y_cross = np.zeros_like(value_axisy)
value_total_y_not   = np.zeros_like(value_axisy)
temp_sum            = np.zeros_like(value_axisy)

value_num   = np.zeros_like(value_axisy)

for i in range(50):
    value_total_x_cross[i] = np.size(workx2d[(value_grid[i]<=gamma[:,-1]) & (value_grid[i+1]>gamma[:,-1]) & (workx2d[:,-1]>worky2d[:,-1]) & (cross_or_not==1) ,-1])
    value_total_x_not[i] = np.size(workx2d[(value_grid[i]<=gamma[:,-1]) & (value_grid[i+1]>gamma[:,-1]) & (workx2d[:,-1]>worky2d[:,-1]) & (cross_or_not==0) ,-1])
    value_total_y_cross[i] = np.size(workx2d[(value_grid[i]<=gamma[:,-1]) & (value_grid[i+1]>gamma[:,-1]) & (workx2d[:,-1]<=worky2d[:,-1]) & (cross_or_not==1) ,-1])
    value_total_y_not[i] = np.size(workx2d[(value_grid[i]<=gamma[:,-1]) & (value_grid[i+1]>gamma[:,-1]) & (workx2d[:,-1]<=worky2d[:,-1]) & (cross_or_not==0) ,-1])
    temp_sum[i] = value_total_x_cross[i]+value_total_x_not[i]+value_total_y_cross[i]+value_total_y_not[i]
    print('x-cross:',value_total_x_cross[i]/temp_sum[i],'x-not:',value_total_x_not[i]/temp_sum[i],'y-cross:',value_total_y_cross[i]/temp_sum[i],'x-not:',value_total_y_not[i]/temp_sum[i])


width=10
pl=plt.bar(value_axisx, value_total_x_cross/temp_sum*100, width, color='lightcoral',edgecolor='black',linewidth=1,label='W$_x$>W$_y$ and w/   crossing')
pl=plt.bar(value_axisx, value_total_x_not/temp_sum*100, width, bottom=value_total_x_cross/temp_sum*100, color='orangered',edgecolor='black',linewidth=1,label='W$_x$>W$_y$ and w/o crossing')
pl=plt.bar(value_axisx, value_total_y_cross/temp_sum*100, width, bottom=(value_total_x_cross+value_total_x_not)/temp_sum*100, color='dodgerblue',edgecolor='black',linewidth=1,label='W$_x$<W$_y$ and w/   crossing')
pl=plt.bar(value_axisx, value_total_y_not/temp_sum*100, width, bottom=(value_total_x_cross+value_total_x_not+value_total_y_cross)/temp_sum*100, color='grey',edgecolor='black',linewidth=1,label='W$_x$<W$_y$ and w/o crossing')

plt.xlim(-10,717)
plt.ylim(0,103)
plt.xlabel('$\epsilon_e$ [m$_e$c$^2$]',fontdict=font)
plt.ylabel('Fraction of LDA and \n TDA electrons',fontdict=font)
plt.xticks(fontsize=22); plt.yticks(fontsize=22);
plt.legend(loc='lower center',fontsize=20,framealpha=0.8,bbox_to_anchor=(0.5, -0.36),ncol=2)
#plt.text(200,650,' t=400fs',fontdict=font)
plt.subplots_adjust(left=0.15, bottom=0.25, right=0.98, top=0.96,
                wspace=None, hspace=None)

#plt.show()
#lt.figure(figsize=(100,100))
fig = plt.gcf()
fig.set_size_inches(12., 8.0)
fig.savefig('./figure_wrap_up/work_cross_x_y'+str(n).zfill(4)+'.png',format='png',dpi=160)
#plt.close("all")

print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
