#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
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

    #### header data ####
    data = sdf.read("./Data_a20/0007.sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y)

    
    line_y1 = np.linspace(-12,-3.2,1001)
    line_x1 = np.zeros_like(line_y1)+5
    
    line_y2 = np.linspace(3.2,12,1001)
    line_x2 = np.zeros_like(line_y2)+5
    
    line_x3 = np.linspace(5,50.0,1001)
    line_y3 = np.zeros_like(line_x3)+3.2
    
    line_x4 = np.linspace(5,50.0,1001)
    line_y4 = np.zeros_like(line_x4)-3.2
    
    
    ax=plt.subplot(3,2,2)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    #den = data['Derived/Number_Density/Subset_high_e'].data/denunit
    
    #levels = np.logspace(-0.3, 2, 51)/2
    #den.T[den.T > 49.999]=49.999 
    #plt.contourf(X, Y, den.T, levels=levels, norm=mcolors.LogNorm(vmin=levels.min(), vmax=levels.max()), cmap='binary')
    ##### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar(ticks=np.logspace(0.0, 2.0, 3))
    #cbar.set_label(r'$n_e\ [n_c]$', fontdict=font)
    
    ex = data['Electric Field/Ey_averaged'].data/exunit
    #ex[ex > 8]=8 
    #ex[ex < -8]=-8 
    print(ex.max())
    #eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 32)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    
    ax.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
    #### manifesting colorbar, changing label and axis properties ####
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -5.0, 0.0, 5.0, 10.0],orientation='horizontal')
    #cbar.set_label('$\overline{E_y}$ [m$_e$c$\omega_0$/e]',fontdict=font)      
    ax.text(21.5,6.6,'$\overline{E}_y$ [m$_e$c$\omega$/|e|]',fontdict=font)      
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.xlabel(r'X [$\lambda_0$]',fontdict=font)
    #ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0)
    ax.tick_params(axis='x',labelsize=0)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)
    
    
    ax=plt.subplot(3,2,4)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    ex = data['Magnetic Field/Bz_averaged'].data/bxunit
    #ex[ex > 8]=8 
    #ex[ex < -8]=-8 
    #eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 32)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    #### manifesting colorbar, changing label and axis properties ####
    ax.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -5.0, 0.0, 5.0, 10.0],orientation='horizontal')
    #cbar.set_label('$\overline{B_z}$ [m$_e\omega_0$/e]',fontdict=font) 
    ax.text(21.5,6.6,'$\overline{B}_z$ [m$_ec\omega$/|e|]',fontdict=font) 
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.xlabel(r'X [$\lambda_0$]',fontdict=font)
    #ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0)
    ax.tick_params(axis='x',labelsize=0)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)
    
    
    ax=plt.subplot(3,2,6)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    ex = data['Magnetic Field/Bz_averaged'].data/bxunit-data['Electric Field/Ey_averaged'].data/exunit
    ex[ex > 5]=5
    ex[ex < -5]=-5 
    #eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 32)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    ax.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -5.0, 0.0, 5.0, 10.0],orientation='horizontal')
   # cbar.set_label('-$\overline{E_y}$+c$\overline{B_z}$ [m$_e$c$\omega_0$/e]',fontdict=font)   
    ax.text(18.5,6.6,'-$\overline{E}_y$+$\overline{B}_z$ [m$_e$c$\omega$/|e|]',fontdict=font)   
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #### manifesting colorbar, changing label and axis properties ####
    #ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
    ax.set_xlabel('X [$\lambda$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=0)
    ax.tick_params(axis='x',labelsize=20)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)


    
    data_e = sdf.read("./Data_uniform/e_fields0006.sdf",dict=True)
    data_b = sdf.read("./Data_uniform/b_fields0006.sdf",dict=True)
    header=data_e['Header']
    time=header['time']
    x  = data_e['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data_e['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y)

    line_y = np.linspace(-12,12,1001)
    line_x = np.zeros_like(line_y)+5
    
    ax=plt.subplot(3,2,1)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    #den = data['Derived/Number_Density/Subset_high_e'].data/denunit
    
    ex = data_e['Electric Field/Ey_averaged'].data/exunit
 #   ex[ex > 8]=8 
 #   ex[ex < -8]=-8 
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 64)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    ax.plot(line_x,line_y,linewidth=3,linestyle=':',color='k')
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -10.0/2, 0.0, 10.0/2, 10.0],orientation='horizontal')
    #cbar.set_label('$\overline{E_y}$ [m$_e$c$\omega_0$/e]',fontdict=font)      
    ax.text(21.5,6.6,'$\overline{E}_y$ [m$_e$c$\omega$/|e|]',fontdict=font)      
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.xlabel(r'X [$\lambda_0$]',fontdict=font)
 #   ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
 #   ax.set_xlabel('X [$\lambda_0$]',fontdict=font)
    ax.set_ylabel('Y [$\lambda$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=20)
    ax.tick_params(axis='x',labelsize=0)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)
    
    ax=plt.subplot(3,2,3)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    ex = data_b['Magnetic Field/Bz_averaged'].data/bxunit
 #   ex[ex > 8]=8 
 #   ex[ex < -8]=-8 
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 64)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    ax.plot(line_x,line_y,linewidth=3,linestyle=':',color='k')
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -10.0/2, 0.0, 10.0/2, 10.0],orientation='horizontal')
    #cbar.set_label('$\overline{B_z}$ [m$_e\omega_0$/e]',fontdict=font)    
    ax.text(21.5,6.6,'$\overline{B}_z$ [m$_ec\omega$/|e|]',fontdict=font)    
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.xlabel(r'X [$\lambda_0$]',fontdict=font)
    #ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
    #ax.set_xlabel('X [$\lambda_0$]',fontdict=font)
    ax.set_ylabel('Y [$\lambda$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=20)
    ax.tick_params(axis='x',labelsize=0)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)
    
    ax=plt.subplot(3,2,5)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper right')
    ex = data_b['Magnetic Field/Bz_averaged'].data/bxunit-data_e['Electric Field/Ey_averaged'].data/exunit
 #   ex[ex > 8]=8 
 #   ex[ex < -8]=-8 
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-10, 10, 64)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-10.0, vmax=10.0), cmap=cm.seismic)
    ax.plot(line_x,line_y,linewidth=3,linestyle=':',color='k')
    cbar=plt.colorbar(image,cax=axin1,ticks=[-10.0, -10.0/2, 0.0, 10.0/2, 10.0],orientation='horizontal')
    #cbar.set_label('-$\overline{E_y}$+c$\overline{B_z}$ [m$_e$c$\omega_0$/e]',fontdict=font)   
    ax.text(18.5,6.6,'-$\overline{E}_y$+$\overline{B}_z$ [m$_e$c$\omega$/|e|]',fontdict=font)   
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #cbar=plt.colorbar()
#    ax.set_ylabel('Y [$\lambda_0$]',fontdict=font)
    ax.set_ylabel('Y [$\lambda$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=20)
    ax.tick_params(axis='x',labelsize=20)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.set_xlabel('X [$\lambda$]',fontdict=font)
    ax.yaxis.set_ticks(ticks=[-8.0,-4.0,0,4.0,8.0])
    ax.text(6,6.6,'t = 117 fs',fontdict=font)
    ax.set_xlim(0,40)
    ax.set_ylim(-8.5,8.5)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.041, hspace=0.051)
    
    fig = plt.gcf()
    fig.set_size_inches(24, 12.8)
    fig.savefig('./figure_wrap_up/Ey_Bz_averaged_1.png',format='png',dpi=160)
    #plt.close("all")
