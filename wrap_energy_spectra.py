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
  wavelength=     1.0e-6
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
      binsize = 200
      en_grid = np.linspace(0.5,799.5,200)
      en_bin  = np.linspace(0,800.0,201)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
      return (en_grid, en_value)

  to_path='./figure_wrap_up/'

  ######### Parameter you should set ###########
  start   =  12  # start time
  stop    =  12  # end time
  step    =  1  # the interval or step
  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  youwant =  ['electron_en']#,'electron_no_en']#,'ey','ex','ey_averaged','bz','bz_averaged','Subset_high_e_density','Subset_high_e_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  if (os.path.isdir('jpg') == False):
    os.mkdir('jpg')
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y)
    
    for name in youwant:
        
        px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
        py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
        grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
        grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
        work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
        work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
        

        gg = (px**2+py**2+1.0)**0.5
        ww = data['Particles/Weight/subset_high_e/electron'].data*24e-6
        theta = np.arctan2(py,px)*180.0/np.pi
    

        px_x_d = px[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
        py_x_d = py[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
        gg_x_d = gg[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
        ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
        dist_x1, den1 = pxpy_to_energy(gg_x_d,ww_x_d)

        px_y_d = px[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
        py_y_d = py[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
        gg_y_d = gg[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
        ww_y_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
        dist_x2, den2 = pxpy_to_energy(gg_y_d,ww_y_d)

        px_d = px[ (grid_x>5) & (abs(grid_y)<3.) ]
        py_d = py[ (grid_x>5) & (abs(grid_y)<3.) ]
        gg_d = gg[ (grid_x>5) & (abs(grid_y)<3.) ]
        ww_d = ww[ (grid_x>5) & (abs(grid_y)<3.) ]
        dist_x, den = pxpy_to_energy(gg_d,ww_d)
        #plt.plot(dist_x*0.51,den,'-k',marker='o',markersize=10,linewidth=3,markevery=8)
        #plt.plot(dist_x1*0.51,den1,'-r',marker='^',markersize=10,linewidth=3,markevery=8)
        #plt.plot(dist_x2*0.51,den2,'-b',marker='s',markersize=10,linewidth=3,markevery=8)
        plt.plot(dist_x*0.51,den,'-k',linewidth=4,label='Total')
        plt.plot(dist_x1*0.51,den1,':r',linewidth=4,label='Work$_x$ > Work$_y$')
        plt.plot(dist_x2*0.51,den2,':b',linewidth=4,label='Work$_x$ < Work$_y$')


        #### manifesting colorbar, changing label and axis properties ####
        plt.xlabel('Energy [MeV]',fontdict=font)
        plt.ylabel('dN/dE [A.U.]',fontdict=font)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);
        plt.yscale('log')
        plt.xlim(0,400)
        plt.legend(loc='best',fontsize=20,framealpha=0.5)
        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
        fig = plt.gcf()
        fig.set_size_inches(10, 7.)
        fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
        plt.close("all")
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  
