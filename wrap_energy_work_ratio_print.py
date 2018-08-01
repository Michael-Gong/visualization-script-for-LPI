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

#choice = np.random.choice(range(px.size), 10000, replace=False)
gamma  = (px**2+py**2+1.)**0.5

ratio=(work_x+work_y+1)/gamma

for i in range(200):
    print('(work_x+work_y+1)/gamma:',(ratio[i]-1)*100.0,'%')


print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
