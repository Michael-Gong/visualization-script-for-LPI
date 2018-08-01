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

directory = './txt_10july/'
xx = np.loadtxt(directory+'xx2d_x.txt')
yy = np.loadtxt(directory+'yy2d_x.txt')
workx2d = np.loadtxt(directory+'workx2d_x.txt')
worky2d = np.loadtxt(directory+'worky2d_x.txt')

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
    



#color_index = abs(theta)

#    plt.subplot()
plt.scatter(workx2d[cross_or_not==0,-1], worky2d[cross_or_not==0,-1], c='crimson', s=20, edgecolors='None', alpha=0.96,label='w/o crossing')
plt.scatter(workx2d[cross_or_not==1,-1], worky2d[cross_or_not==1,-1], c='royalblue', s=20, edgecolors='None', alpha=0.96,label='w/   crossing')
#plt.scatter(workx2d[cross_or_not==1,-1], worky2d[cross_or_not==1,-1], c='royalblue', s=17, edgecolors='black', alpha=0.95,label='cross')
#plt.scatter(workx2d[cross_or_not==0,-1], worky2d[cross_or_not==0,-1], c='crimson', s=17, edgecolors='black', alpha=0.95,label='not_cross')
plt.legend(loc='upper right',scatterpoints=1,markerscale=3, framealpha=0.5,fontsize=27)

plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
plt.xlim(-250,750)
plt.ylim(-250,750)
plt.xlabel('W$_x$ [m$_e$c$^2$]',fontdict=font)
plt.ylabel('W$_y$ [m$_e$c$^2$]',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
plt.text(-100,650,' t = 400 fs',fontdict=font)
plt.subplots_adjust(left=0.18, bottom=None, right=0.97, top=None,
                wspace=None, hspace=None)

#plt.show()
#lt.figure(figsize=(100,100))
fig = plt.gcf()
fig.set_size_inches(10.5, 11)
fig.savefig('./figure_wrap_up/cross_or_not.png',format='png',dpi=160)
#plt.close("all")

print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
