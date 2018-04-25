import sdf
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
#from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
#from optparse import OptionParser
#import os
#from colour import Color

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
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  
######### Parameter you should set ###########
start   =  0  # start time
stop    =  19  # end time
step    =  1  # the interval or step

#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']

    px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
    py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
    grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
    grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
    work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
    work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
    gg = (px**2+py**2+1)**0.5
    
    px = px[gg > 2.0]
    py = py[gg > 2.0]
    grid_x = grid_x[gg > 2.0]
    grid_y = grid_y[gg > 2.0]
    work_x = work_x[gg > 2.0]
    work_y = work_y[gg > 2.0]

    np.savetxt('./txt/px_'+str(n).zfill(4)+'sdf.txt',px)
    np.savetxt('./txt/py_'+str(n).zfill(4)+'sdf.txt',py)
    np.savetxt('./txt/grid_x_'+str(n).zfill(4)+'sdf.txt',grid_x)
    np.savetxt('./txt/grid_y_'+str(n).zfill(4)+'sdf.txt',grid_y)
    np.savetxt('./txt/work_x_'+str(n).zfill(4)+'sdf.txt',work_x)
    np.savetxt('./txt/work_y_'+str(n).zfill(4)+'sdf.txt',work_y)

    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')