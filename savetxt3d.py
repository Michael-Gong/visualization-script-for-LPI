import sdf 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal



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
print 'electric field unit: '+str(exunit)
print 'magnetic field unit: '+str(bxunit)
print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  



######### Parameter you should set ###########
start   =  4  # start time
stop    =  4  # end time
step    =  1  # the interval or step

#youwant = ['']
youwant =  ['ey','electron_density']
#youwant.append('carbon_ekbar')
#youwant.append('positron_ekbar')
#youwant.append('electron_en')
#youwant.append('photon_en')
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...
#youwant dist_fn electron_x_px...

######### Script code drawing figure ################
for n in range(start,stop+step,step):
  #### header data ####
  data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  X, Y = np.meshgrid(x, y)
  np.savetxt('./txt/x.txt', x)
  np.savetxt('./txt/y.txt', y)
  for name in youwant:
    if (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
              eexx = data['Electric Field/'+str.capitalize(name)].data/exunit
              n3d=len(eexx[0,0,:])
              ex = (eexx[:,:,n3d/2-1]+eexx[:,:,n3d/2])/2 
              np.savetxt('./txt/'+name[0:2]+'.txt', ex)
    elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
              eexx = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
              n3d=len(eexx[0,0,:])
              ex = (eexx[:,:,n3d/2-1]+eexx[:,:,n3d/2])/2
              np.savetxt('./txt/'+name[0:2]+'.txt', ex)
    elif (name[-7:] == 'density'):
              ddeen = data['Derived/Number_Density/'+name[0:-8]].data/denunit
              n3d=len(ddeen[0,0,:])
              den = (ddeen[:,:,n3d/2-1]+ddeen[:,:,n3d/2])/2
              np.savetxt('./txt/density_'+name[0:-8]+'.txt', den)
    elif (name[-5:] == 'ekbar'):
              ddeen = data['Derived/EkBar/'+name[0:-6]].data/(q0*1.0e6)
              n3d=len(ddeen[0,0,:])
              den = (ddeen[:,:,n3d/2-1]+ddeen[:,:,n3d/2])/2
              np.savetxt('./txt/ekbar_'+name[0:-6]+'.txt', den)
  print 'finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%'


