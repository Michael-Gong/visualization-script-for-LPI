#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from matplotlib.colors import ListedColormap

  
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
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
  
  font = {'family' : 'helvetica',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
  
  
  
  n=47
  data = sdf.read("./"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  y  = y[600:1799]
  X, Y = np.meshgrid(x, y)

  plt.subplot(2,1,1)
  ex = data['Electric Field/Ex'].data/exunit
  ex = ex[:,600:1799]
  eee=np.max([-np.min(ex.T),np.max(ex.T)])
  levels = np.linspace(-eee, eee, 40)
  plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdBu)
  #### manifesting colorbar, changing label and axis properties ####
  cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee],orientation='horizontal')
  cbar.set_label(r'$E_x$ [$m_ec\omega/e$]',fontdict=font)        
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('Y [$\mu m$]',fontdict=font)
  plt.xticks(fontsize=20); plt.yticks(fontsize=20);
  plt.title('Ex at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                
  plt.subplot(2,1,2)
  ex = data['Electric Field/Ey'].data/exunit
  ex = ex[:,600:1799]
  eee=np.max([-np.min(ex.T),np.max(ex.T)])
  levels = np.linspace(-eee, eee, 40)
  plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdBu, alpha=0.5)
  cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee],orientation='horizontal')
  cbar.set_label(r'$E_y$ [$m_ec\omega/e$]', fontdict=font)

##generate the transparent colorbar
  cmap = plt.cm.Greys
  my_cmap = cmap(np.arange(cmap.N))
  my_cmap[:,-1] = np.sqrt(np.linspace(0.0, 1, cmap.N))
  my_cmap = ListedColormap(my_cmap)

  den = data['Derived/Number_Density/electron'].data/denunit
  den = den[:,600:1799]
  levels = np.logspace(-3, -1, 40) 
  plt.contourf(X, Y, den.T, levels=levels, cmap=my_cmap)
  #### manifesting colorbar, changing label and axis properties ####
  cbar=plt.colorbar(ticks=np.logspace(-3, -1, 5),orientation='horizontal')
  cbar.set_label(r'$n_e$ [$n_c$]', fontdict=font)
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('Y [$\mu m$]',fontdict=font)
  plt.xticks(fontsize=20); plt.yticks(fontsize=20);
  plt.title('Density at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)

  fig = plt.gcf()
  fig.set_size_inches(12, 24)
  fig.savefig('./field'+str(n).zfill(4)+'.png',format='png',dpi=160)
  plt.close("all")
