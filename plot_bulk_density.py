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
          'size'   : 20,  
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

    for n in range(16,17):
        data_0 = sdf.read("./Data_a20_130_1/"+str(1).zfill(4)+".sdf",dict=True)
        density_e0_0 = data_0['Derived/Number_Density_averaged/electron'].data
        density_e1_0 = data_0['Derived/Number_Density_averaged/electron_no'].data
        density_c0_0 = data_0['Derived/Number_Density_averaged/carbon'].data



        data = sdf.read("./Data_a20_130_1/"+str(n).zfill(4)+".sdf",dict=True)
        x    = data['Grid/Grid_mid'].data[0]/1e-6
        y    = data['Grid/Grid_mid'].data[1]/1e-6
        X, Y = np.meshgrid(x,y)
    
        density_e0 = data['Derived/Number_Density_averaged/electron'].data-density_e0_0
        density_e1 = data['Derived/Number_Density_averaged/electron_no'].data-density_e1_0
        density_c0 = data['Derived/Number_Density_averaged/carbon'].data-density_c0_0
     
        charge_density = (density_e0+density_e1)*(-1.0)+density_c0*6.0
    
    #    charge_density = np.sum(charge_density[:,320:386] ,axis=1)/66/denunit
        
        
        plt.subplot(1,1,1)
        #ex = data['Electric Field/Ex_averaged'].data/exunit
        #ex=ex-ex_ave
     #   ex[ex > 8]=8 
     #   ex[ex < -8]=-8 
        den = charge_density/denunit
        eee=np.max([-np.min(den),np.max(den)])
        half_a = -np.linspace(0.01,2.0,32)
        half_b = np.linspace(0,30,32)
        levels = np.concatenate((half_a[::-1],half_b),axis=0)
        #plt.contourf(X, Y, den.T, levels=levels, norm=colors.SymLogNorm(linthresh=1, linscale=0.1, vmin=-1e1, vmax=1e2, midpoint=0), cmap=cm.bwr)
        plt.contourf(X, Y, den.T, levels=levels, norm=MidpointNormalize(midpoint=0.), cmap=cm.bwr)
        print('maximum:',eee)
    #    line_x  =  np.linspace(5,200,1001)
    #    line_y1 =  np.zeros_like(line_x)+3.2
    #    line_y2 =  np.zeros_like(line_x)-3.2
    #    ax.plot(line_x,line_y1,linewidth=3,linestyle=':',color='k')
    #    ax.plot(line_x,line_y2,linewidth=3,linestyle=':',color='k')
        #plt.plot(line_x,line_y,linewidth=3,linestyle=':',color='k')
        #plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.seismic)
        #### manifesting colorbar, changing label and axis properties ####
        #cbar=plt.colorbar()
        cbar=plt.colorbar(ticks=[ -2, -1, 0.0, 15, 30])
        cbar.set_label('$n^-$+$n^+$ [$n_c$]',fontdict=font2)
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)
        #plt.text(30,7,r'$E_x\ [m_ec\omega_0/e]$',fontdict=font)
        plt.xlabel(r'X [$\lambda$]',fontdict=font)
        plt.ylabel(r'Y [$\lambda$]',fontdict=font)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(0,100)
        plt.ylim(-4,4)
    #    ax.text(60,8,'t = '+'400 fs',fontdict=font)
    #    ax.text(90,8,'maximum = '+str(round(eee,2)),fontdict=font)
        #plt.xticks(fontsize=20); plt.yticks([-10,-5,0,5,10],fontsize=20);
        #plt.xticks([])
        #plt.title('At '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
        #ax.set_xlim(87+(n-22)*5,123+(n-22)*5)
    #    ax.set_xlim(0,120)
    #    ax.set_ylim(-12,12)
        #plt.yticks(np.linspace(-8,8,5))
    
        plt.subplots_adjust(left=None, bottom=0.25, right=None, top=None,
                    wspace=0.04, hspace=None)
    
        fig = plt.gcf()
        fig.set_size_inches(15, 3.)
        fig.savefig('./jpg_a20_130_1/caculate_density_'+str(n).zfill(4)+'.png',format='png',dpi=160)
        plt.close("all")
        print('finished ',n)

        grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
        grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
        px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
        py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
        gg = (px**2+py**2)**0.5
        grid_x = grid_x [(abs(grid_y) < 3.2) & (gg > 2.0)]
        px = px [(abs(grid_y) < 3.2) & (gg > 2.0)]
        py = py [(abs(grid_y) < 3.2) & (gg > 2.0)]
        gg = (px**2+py**2+1)**0.5
        plt.subplot(1,1,1)
        plt.scatter(grid_x, gg, c='blue', s=5, edgecolors='None', alpha=1.)

        plt.xlim(0,100)
        #    plt.ylim(0,400)
        #plt.xlabel(r'X [$\lambda$]',fontdict=font)
        plt.ylabel(r'$\gamma$',fontdict=font)
        plt.xticks(fontsize=0); 
        plt.yticks([0,300,600],fontsize=20);
        plt.ylim(0,600.0)
        plt.subplots_adjust(left=None, bottom=0.25, right=None, top=None,
                    wspace=0.04, hspace=None)
    
        fig = plt.gcf()
        fig.set_size_inches(12, 2.8)
        fig.savefig('./jpg_a20_130_1/caculate_scatter_'+str(n).zfill(4)+'.png',format='png',dpi=160)
        plt.close("all")
        print('finished ',n)


 
