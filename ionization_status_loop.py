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
      binsize = 500
      en_grid = np.linspace(1,999,500)
      en_bin  = np.linspace(0,1000.0,501)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(1000.0/binsize)
      return (en_grid, en_value)
  
  #from_list = ['./laser_a001/','./laser_a002/','./laser_a005/','./laser_a010/','./laser_a020/','./laser_a050/','./laser_a100/','./laser_a200/','./laser_a500/']
#  from_list = ['./laser_a1000/','./laser_a2000/','./laser_a5000/']
  #from_list = ['./a1/','./a2/','./a3/','./a4/']
  from_list = ['./carbon_ion_a05/']

  for from_path in from_list:
    #  from_path='./a4/'
      to_path = './jpg_carbon_ion_a05/'
      ######### Parameter you should set ###########
      start   =  0  # start time
      stop    =  24  # end time
      step    =  1  # the interval or step
      
    #  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
      youwant =  ['C','C1','C2','C3','C4','C5','C6']
      #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
      #youwant Derived electron_density,electron_ekbar...
      #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
      #if (os.path.isdir('jpg') == False):
      #  os.mkdir('jpg')
      ######### Script code drawing figure ################
      num_ion = np.zeros([np.size(youwant),stop-start+step])
      frac_ion= np.zeros_like(num_ion)
      time_grid = np.zeros_like(num_ion[0,:])
      for n in range(start,stop+step,step):
          data = sdf.read(from_path+'tot'+str(n).zfill(4)+'.sdf',dict=True)
          header=data['Header']
          time=header['time']
          time_grid[n]=time
          for i in range(np.size(youwant)):
              if 'Particles/ID/subset_sub/'+youwant[i] in data:
                  num_ion[i,n-start] = np.size(data['Particles/ID/subset_sub/'+youwant[i]].data)
          print('finish '+str(n).zfill(4))    
          frac_ion[:,n-start] = num_ion[:,n-start]/np.sum(num_ion[:,n-start])    
    
      for i in range(np.size(youwant)):
          plt.plot(time_grid/1e-15,100.0*frac_ion[i,:],linewidth=2.5, label=youwant[i])
        
      #### manifesting colorbar, changing label and axis properties ####
      plt.xlabel('time [fs]',fontdict=font)
      plt.ylabel('fraction [%]',fontdict=font)
      plt.xticks(fontsize=20); plt.yticks(fontsize=20);
      #plt.yscale('log')
      plt.ylim(0,50)
      #plt.xlim(5,600)
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
    
      plt.legend(loc='best',fontsize=10,framealpha=0.0)
      plt.subplots_adjust(left=None, bottom=0.15, right=0.95, top=0.95,
                    wspace=None, hspace=None)
        #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(10.0, 8.5)
      fig.savefig(to_path+'ion_fraction.png',format='png',dpi=160)
      plt.close("all")
