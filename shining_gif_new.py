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
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
  
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
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('salmon')
  c_blue= matplotlib.colors.colorConverter.to_rgba('skyblue')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128) 
##end for transparent colorbar##
 
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

  ######### Parameter you should set ###########
  start   =  300  # start time
  stop    =  499  # end time
  step    =  1  # the interval or step
  from_path='../conductor/Data_new/'
  from_path1='./jpg_shining_g400/'
  to_path  ='./jpg_shining_new/'

  xx_2d_x = np.loadtxt(from_path1+'xx_2d_x.txt') 
  yy_2d_x = np.loadtxt(from_path1+'yy_2d_x.txt') 
  gg_2d_x = np.loadtxt(from_path1+'gg_2d_x.txt') 

  xx_2d_y = np.loadtxt(from_path1+'xx_2d_y.txt') 
  yy_2d_y = np.loadtxt(from_path1+'yy_2d_y.txt') 
  gg_2d_y = np.loadtxt(from_path1+'gg_2d_y.txt') 
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    name = 'Subset_high_e_density'
    den = data['Derived/Number_Density/electron'].data/denunit
    den_p = data['Derived/Number_Density/electron_no'].data/denunit
    den = den+den_p
    den = den[:,:]
    
    if np.min(den.T) == np.max(den.T):
            continue
    levels = np.linspace(0.0, 52.499, 101)
    den.T[den.T > 52.499]=52.499 
    plt.contourf(X, Y, den.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap='Greys')
    #### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar(ticks=np.linspace(0.0, 50.0, 5))
    #cbar.set_label(name+'[$n_c$]', fontdict=font)

    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex[:,:]
    if np.min(ex.T) == np.max(ex.T):
            continue
    levels = np.linspace(-22.5, 22.5, 11)
    ex.T[ex.T < -22.499]=-22.499
    ex.T[ex.T >  22.499]= 22.499
    plt.contour(X, Y, ex.T, levels=levels, cmap=cmap_br)
    #### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar(ticks=np.linspace(-22.5, 22.5, 5))
    #cbar.set_label('Normalized electric field',fontdict=font)        
    if n-50 >= 6:
        print('n-50=',n-50)
        index=np.arange(n-50-6,n-50,1)
        print('index is',index)
        select = np.where(gg_2d_x[:,n-50] > 1.0)
        print('select is',select)
        print('xx_2d_x[:,index][select,:] shape is',xx_2d_x[:,index][select,:].shape)
        #print('yy_2d_x[select,index] shape is',(yy_2d_x[select,index]).shape)
        #print('xx_2d_x[select,0].size is',xx_2d_x[select,0].size)
        plt.scatter(xx_2d_x[:,index][select,:].reshape(np.array(select).size,index.size), yy_2d_x[:,index][select,:].reshape(np.array(select).size,index.size), c=np.tile(np.arange(index.size)*5, (xx_2d_x[select,0].size, 1)), s=np.tile(np.arange(index.size)*5, (xx_2d_x[select,0].size, 1)), cmap='rainbow', edgecolors='None')
        select = np.where(gg_2d_y[:,n-50] > 1.0)
        print('select is',select)
        print('xx_2d_y[:,index][select,:] shape is',xx_2d_y[:,index][select,:].shape)
        #print('yy_2d_y[select,index] shape is',(yy_2d_y[select,index]).shape)
        #print('xx_2d_y[select,0].size is',xx_2d_y[select,0].size)
        plt.scatter(xx_2d_y[:,index][select,:].reshape(np.array(select).size,index.size), yy_2d_y[:,index][select,:].reshape(np.array(select).size,index.size), c=np.tile(np.arange(index.size)*5, (xx_2d_y[select,0].size, 1)), s=np.tile(np.arange(index.size)*5, (xx_2d_y[select,0].size, 1)), cmap='winter_r', edgecolors='None')
    plt.ylim(-5,5)
    plt.xlim(0+(n-299.0)*0.1,30+(n-299.0)*0.1)
    plt.xlabel('X [$\lambda_0$]',fontdict=font)
    plt.ylabel('Y [$\lambda_0$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);

    fig = plt.gcf()
    fig.set_size_inches(24, 7.0)
    fig.savefig(to_path+str(n).zfill(4)+'.png',format='png',dpi=50)
    plt.close("all")


    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

