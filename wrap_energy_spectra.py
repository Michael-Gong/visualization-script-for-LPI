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
  
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     1.06e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
  
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####


  def pxpy_to_energy(gamma, weight):
      binsize = 600
      en_grid = np.linspace(0.5,599.5,600)
      en_bin  = np.linspace(0,600.0,601)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
      return (en_grid, en_value)

  to_path='./'

  ######### Parameter you should set ###########
  start   =  12  # start time
  stop    =  30  # end time
  step    =  1  # the interval or step
  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  youwant =  ['electron_en']#,'electron_no_en']#,'ey','ex','ey_averaged','bz','bz_averaged','Subset_high_e_density','Subset_high_e_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  weight = V*denunit*n0/50.0  
  
  data = sdf.read("./i_tot_loc"+str(17).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  gg = (px**2+py**2+py**2+1)**0.5
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi

  ek_1 = (gg - 1.0)*0.51*1836
  ww_1 = np.zeros_like(ek_1) + weight
  dist_x1, den1 = pxpy_to_energy(ek_1,ww_1)

  ek_2 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
  ww_2 = np.zeros_like(ek_2) + weight
  dist_x2, den2 = pxpy_to_energy(ek_2,ww_2)


  data = sdf.read("./i_tot_loc"+str(27).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  gg = (px**2+py**2+py**2+1)**0.5
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi

  ek_3 = (gg - 1.0)*0.51*1836
  ww_3 = np.zeros_like(ek_3) + weight
  dist_x3, den3 = pxpy_to_energy(ek_3,ww_3)

  ek_4 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
  ww_4 = np.zeros_like(ek_4) + weight
  dist_x4, den4 = pxpy_to_energy(ek_4,ww_4)


  plt.plot(dist_x1,den1,':b',linewidth=4, label='total')
  plt.plot(dist_x2,den2,':r',linewidth=4, label=r'$\theta$'+'< 10$^o$')
  plt.plot(dist_x3,den3,'-b',linewidth=4, label='total')
  plt.plot(dist_x4,den4,'-r',linewidth=4, label=r'$\theta$'+'< 10$^o$')

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('Energy [MeV]',fontdict=font)
  plt.ylabel('dN/dE [per MeV]',fontdict=font)
  plt.xticks(fontsize=20); plt.yticks(fontsize=20);
  plt.yscale('log')
  plt.ylim(1e7,1e10)
  plt.legend(loc='best',fontsize=20,framealpha=0.5)
#        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(8.0, 6.5)
  fig.savefig('./maybe.png',format='png',dpi=160)
  plt.close("all")
  
