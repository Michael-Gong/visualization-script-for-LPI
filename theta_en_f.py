#!/public/home/users/bio001/tools/python-2.7.11/bin/python
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
start   =  1  # start time
stop    =  1  # end time
step    =  1  # the interval or step

youwant = ['electron_theta_en']#['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
#youwant =  ['bz','ex','ey_averaged','ez','electron_density','carbon_density','photon_density','positron_density','electron_ekbar','photon_ekbar','electron_x_px']
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...
#youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...

######### Script code drawing figure ################
for n in range(start,stop+step,step):
  #### header data ####
  data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
#  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
#  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
#  X, Y = np.meshgrid(x, y)
  
  for name in youwant:
    if (name[-8:] == 'theta_en'):
              denden = data['dist_fn/theta_en/'+name[0:-9]].data[:,:,0]
              den = np.log(denden+1.0)
              if np.min(den.T) == np.max(den.T):
                  continue
              min_value=np.min(den.T[den.T > 0])
              max_value=np.max(den.T)
              levels = np.linspace(min_value, max_value, 40)
              dist_x  = data['Grid/theta_en/'+name[0:-9]].data[0]/np.pi*180
              dist_y  = data['Grid/theta_en/'+name[0:-9]].data[1]/(q0*1.0e6)
              dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
              plt.subplot(2,1,1)
              plt.contourf(dist_X, dist_Y, den.T, levels=levels, cmap=cm.nipy_spectral)
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(ticks=np.linspace(min_value, max_value, 5))
              cbar.set_label(r'$log_{10}\frac{dN}{d\theta dE}$ [A.U.]', fontdict=font)
              plt.xlabel(r'$\theta$ [degree]',fontdict=font)
              plt.ylabel('Energy [MeV]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              xmin, xmax = plt.xlim()
              plt.subplot(2,1,2)
              #plt1 = plt.twinx()
              plt.plot(dist_x,np.sum(denden,axis=1),'-r',linewidth=2.5)
              #plt.xlim(-180,180)
              #plt1.set_ylabel('Normalized '+name)  
              cbar=plt.colorbar(ticks=np.linspace(min_value, max_value, 5))
              plt.xlabel(r'$\theta$ [degree]',fontdict=font)
              plt.ylabel(r'$\frac{dN}{d\theta}$ [A.U.]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.xlim(xmin,xmax)

              fig = plt.gcf()
              fig.set_size_inches(12, 14)
              fig.savefig('./jpg/'+name+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
  print 'finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%'

