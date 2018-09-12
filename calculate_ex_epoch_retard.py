#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec


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
    font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 15,  
          }  
    
##below is for generating mid transparent colorbar
    c_red = matplotlib.colors.colorConverter.to_rgba('red')
    c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
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
    nx = 51
    ny = 36

    ex_field       = np.zeros([36,51])
    ey_field       = np.zeros([36,51])
    grid_x         = np.linspace(90,95,51)*1.e-6
    grid_y         = np.linspace(0,3.5,36)*1.e-6

    data = sdf.read("../Data_a20_130/0022.sdf",dict=True)
    x    = data['Grid/Grid_mid'].data[0]
    y    = data['Grid/Grid_mid'].data[1]
    Grid_x,Grid_y  = np.meshgrid(x,y)

    density_e0 = data['Derived/Number_Density/electron'].data
    density_e1 = data['Derived/Number_Density/electron_no'].data
    density_c0 = data['Derived/Number_Density/carbon'].data

    g0_e       = data['Derived/EkBar/electron'].data/(m0*v0**2)+1.0
    g0_c       = data['Derived/EkBar/carbon'].data/(6*1836*m0*v0**2)+1.0

    for ix in range(nx):
        for iy in range(ny):
            Rx = grid_x[ix] - Grid_x + 1.0e-50
            Ry = grid_y[iy] - Grid_y + 1.0e-50
            g0 = g0_e.T
            vv = (1.0-1/g0**2)**0.5*v0
            t0 = (2*Rx*vv+(4*Rx**2*vv**2+4*(Rx**2+Ry**2)*(v0**2-vv**2))**0.5)/2./(v0**2-vv**2)
            charge_density = (density_e0+density_e1)*(-1.0)
            ex_temp = q0*charge_density.T/2/pi/epsilon0*((grid_x[ix] - Grid_x)/v0/t0)/g0**2/(1.0-(vv**2*t0**2+v0**2*t0**2-Rx**2-Ry**2)/2/v0**2/t0**2)**3/v0/t0
            ey_temp = q0*charge_density.T/2/pi/epsilon0*((grid_y[iy] - Grid_y)/v0/t0)/g0**2/(1.0-(vv**2*t0**2+v0**2*t0**2-Rx**2-Ry**2)/2/v0**2/t0**2)**3/v0/t0
            g0 = g0_c.T
            vv = (1.0-1/g0**2)**0.5*v0
            t0 = (2*Rx*vv+(4*Rx**2*vv**2+4*(Rx**2+Ry**2)*(v0**2-vv**2))**0.5)/2./(v0**2-vv**2)
            charge_density = (density_c0)*6.0
            ex_temp = ex_temp + q0*charge_density.T/2/pi/epsilon0*((grid_x[ix] - Grid_x)/v0/t0)/g0**2/(1.0-(vv**2*t0**2+v0**2*t0**2-Rx**2-Ry**2)/2/v0**2/t0**2)**3/v0/t0
            ey_temp = ey_temp + q0*charge_density.T/2/pi/epsilon0*((grid_y[iy] - Grid_y)/v0/t0)/g0**2/(1.0-(vv**2*t0**2+v0**2*t0**2-Rx**2-Ry**2)/2/v0**2/t0**2)**3/v0/t0

            ex_field[iy,ix] = sum(sum(ex_temp))*(x[-1]-x[-2])*(y[-1]-y[-2])
            ey_field[iy,ix] = sum(sum(ey_temp))*(x[-1]-x[-2])*(y[-1]-y[-2])

        print('total ix is ',nx,'; we finish:',ix)
   
    np.savetxt('./ex_field_retard.txt',ex_field)
    np.savetxt('./ey_field_retard.txt',ey_field)
