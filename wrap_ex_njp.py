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

    #### header data ####
    data = sdf.read("./0012_laser.sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y)
    
    ax=plt.subplot(1,2,1)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper left')
    ex = data['Electric Field/Ex'].data/exunit
 #   ex[ex > 8]=8 
 #   ex[ex < -8]=-8 
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-2, 2, 64)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-2.0, vmax=2.0), cmap=cm.jet)
    line_x  =  np.linspace(5,200,1001)
    R_length = 3.14*4**2/1.0
    line_y1 = 4.0*(1+((line_x-5.0)/R_length)**2.0)**0.5
    line_y2 =-4.0*(1+((line_x-5.0)/R_length)**2.0)**0.5
    ax.plot(line_x,line_y1,linewidth=3,linestyle=':',color='grey')
    ax.plot(line_x,line_y2,linewidth=3,linestyle=':',color='grey')
    #plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.seismic)
    #### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar()
    cbar=plt.colorbar(image,cax=axin1,ticks=[-2.0, -2.0/2, 0.0, 2.0/2, 2.0],orientation='horizontal')
    cbar.set_label('E$_x$ [m$_e$c$\omega$/|e|]',fontdict=font2)    
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.text(30,7,r'$E_x\ [m_ec\omega_0/e]$',fontdict=font)
    ax.set_xlabel(r'X [$\lambda$]',fontdict=font)
    ax.set_ylabel(r'Y [$\lambda$]',fontdict=font)
    ax.tick_params(axis='y',labelsize=20)
    ax.tick_params(axis='x',labelsize=20)
    #plt.xticks([])
    ax.text(110,17.5,'t = 400 fs',fontdict=font)
    ax.set_xlim(87,123)
    ax.set_ylim(-22.5,22.5)
    #plt.yticks(np.linspace(-8,8,5))
    
    
    data = sdf.read("./Data_a20_130/0022.sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6+10.0
    print('ok')
    y  = np.linspace(-31.975,31.975,1280)
    X, Y = np.meshgrid(x, y)
    
    assist = np.zeros([12500,1280])
    ax=plt.subplot(1,2,2)
    axin1 = inset_axes(ax, width='20%', height='5%', loc='upper left')
    et = data['Electric Field/Ex'].data/exunit
    assist[:,320:960]=et
    tempt1 = et[:,-320:]
    tempt2 = et[:,:320]
    assist[:,0:320]=tempt1+assist[:,0:320]
    assist[:,960:]=tempt2+assist[:,960:]
    ex = assist 
    print('ex shape is',ex.shape)
    #ex = data['Electric Field/Ex_averaged'].data/exunit
    #ex=ex-ex_ave
 #   ex[ex > 8]=8 
 #   ex[ex < -8]=-8 
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-2, 2, 64)
    image=ax.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-2.0, vmax=2.0), cmap=cm.jet)
    line_x  =  np.linspace(5,200,1001)
    line_y1 =  np.zeros_like(line_x)+3.2
    line_y2 =  np.zeros_like(line_x)-3.2
    ax.plot(line_x,line_y1,linewidth=3,linestyle=':',color='k')
    ax.plot(line_x,line_y2,linewidth=3,linestyle=':',color='k')
    #plt.plot(line_x,line_y,linewidth=3,linestyle=':',color='k')
    #plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.seismic)
    #### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar()
    cbar=plt.colorbar(image,cax=axin1,ticks=[-2.0, -2.0/2, 0.0, 2.0/2, 2.0],orientation='horizontal')
    cbar.set_label('E$_x$ [m$_e$c$\omega$/|e|]',fontdict=font2)    
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=15)
    #plt.text(30,7,r'$E_x\ [m_ec\omega_0/e]$',fontdict=font)
    ax.set_xlabel(r'X [$\lambda$]',fontdict=font)
    #ax.set_ylabel(r'Y [$\lambda_0$]',fontdict=font)
    ax.tick_params(axis='x',labelsize=20)
    ax.tick_params(axis='y',labelsize=0)
    ax.text(110,17.5,'t = '+'400 fs',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks([-10,-5,0,5,10],fontsize=20);
    #plt.xticks([])
    #plt.title('At '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
    ax.set_xlim(87,123)
    ax.set_ylim(-22.5,22.5)
    #plt.yticks(np.linspace(-8,8,5))
    
    
    
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.04, hspace=None)
    
    fig = plt.gcf()
    fig.set_size_inches(24, 6.4)
    fig.savefig('./figure_wrap_up/Ex_1_njp_1.png',format='png',dpi=160)
    plt.close("all")
