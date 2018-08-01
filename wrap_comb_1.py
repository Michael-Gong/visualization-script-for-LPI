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
          'size'   : 25,  
          }  
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 12,  
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

 
  
  ######### Parameter you should set ###########
  start   =  210  # start time
  stop    =  210  # end time
  step    =  1  # the interval or step
  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read("../conductor/Data_new/"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    x  = x[0:2400]
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    name = 'Subset_high_e_density'
    den = data['Derived/Number_Density/electron'].data/denunit
    den_p = data['Derived/Number_Density/electron_no'].data/denunit
    den = den+den_p
    den = den[0:2400,:]
    
    if np.min(den.T) == np.max(den.T):
            continue
    levels = np.linspace(0.0, 52.499, 101)
    den.T[den.T > 52.499]=52.499 

    #gs = gridspec.GridSpec(2, 2, width_ratios=[6, 1], height_ratios=[1, 3])

    
    ax=plt.subplot(1,2,1)
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left',pad=0.2)
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='lower left',pad=0.2)
    #axin1 = inset_axes(ax,width="5%",height="45%",loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #axin2 = inset_axes(ax,width="5%",height="45%",loc='lower right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
    #### manifesting colorbar, changing label and axis properties ####
    image1=ax.contourf(X, Y, den.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap='Greys')
    #cbar=plt.colorbar(image1,cax=axin1,ticks=np.linspace(0.0, 50.0, 5),orientation="horizontal")
    #cbar.set_label('$n_e$ [$n_c$]', fontdict=font2)

    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex[0:2400,:]
    if np.min(ex.T) == np.max(ex.T):
            continue
    levels = np.linspace(-22.5, 22.5, 50)
    ex.T[ex.T < -22.499]=-22.499
    ex.T[ex.T >  22.499]= 22.499
    image2=ax.contourf(X, Y, ex.T, levels=levels, cmap=cmap_br)
    #### manifesting colorbar, changing label and axis properties ####
    #cbar=plt.colorbar(image2,cax=axin2,ticks=np.linspace(-20, 20, 5),orientation="horizontal")
    #cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font2)        
    #ax.text(21.,1.75,'t = '+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    ax.text(21.,1.75,'t = 70 fs',fontdict=font)
    if 'Particles/Px/subset_high_e/electron' in data:
    
        px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
        py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
        grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
        grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
        gg = (px**2+py**2+1)**0.5
        grid_x = grid_x[gg > 1.5]
        grid_y = grid_y[gg > 1.5]
        gg = gg[gg > 1.5] 
        if gg.size > 1.5:
            ppp_i = np.random.choice(gg.size,gg.size,replace=False)
            ppp_x = grid_x[ppp_i]
            ppp_y = grid_y[ppp_i]
            ppp_px = px[ppp_i]
            ppp_py = py[ppp_i]
            ppp_g = gg[ppp_i]
            #plt.scatter(ppp_x, ppp_y, c=ppp_g, s=20, cmap='summer', edgecolors='None', alpha=0.1)
            ax.scatter(ppp_x[abs(ppp_y)<=3.2], ppp_y[abs(ppp_y)<=3.2], s=15, c='green', norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), edgecolors='None', alpha=0.5)
            #plt.scatter(ppp_x[abs(ppp_y)>3.2], ppp_y[abs(ppp_y)>3.2], s=4, c='cyan', norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), edgecolors='None', alpha=0.5)
            #cbar=plt.colorbar( ticks=np.linspace(np.min(ppp_g), np.max(ppp_g), 5) )
            #cbar.set_label(r'$\gamma$',fontdict=font)
    
    ax.set_ylim(-6.,6.)
    ax.set_xlim(0,26)
    ax.set_xlabel('X [$\lambda$]',fontdict=font)
    ax.set_ylabel('Y [$\lambda$]',fontdict=font)
    ax.tick_params(axis='both',labelsize=25) 
    #ax.set_xticklabels(xticklabels,fontdict=font)
    #ax.set_yticklabels(yticklabels,fontdict=font)

#    plt.subplot(gs[0])
#    plt.scatter(ppp_x[abs(ppp_y)<=3.2],ppp_px[abs(ppp_y)<=3.2],s=10,c='green',edgecolors='None',alpha=0.5)
#    plt.xlim(0,30)
#    plt.ylabel('p$_x$ [m$_e$c]', fontdict=font)
#    plt.xticks([])
#    plt.yticks(fontsize=20)
#  
#    plt.subplot(gs[3])
#    plt.scatter(ppp_py[abs(ppp_y)<=3.2],ppp_y[abs(ppp_y)<=3.2],s=10,c='green',edgecolors='None',alpha=0.5)
#    plt.ylim(-6.5,6.5)
#    plt.xlabel('p$_y$ [m$_e$c]', fontdict=font)
#    plt.yticks([])
#    plt.xticks(fontsize=20)
#
#    
#    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.011, hspace=0.051)
#

    #fig=plt.subplot(gs[1])
    #ax1 = fig.add_axes([0.05, 0.85, 0.9, 0.10])
    #ax2 = fig.add_axes([0.05, 0.35, 0.9, 0.10])
    #cmap = mpl.cm.rainbow
    #norm = mpl.colors.Normalize(vmin=0.0, vmax=50)
    #cb1 = mpl.colorbar.ColorbarBase(ax1, cmap='Greys',
    #                            norm=norm,
    #                            orientation='horizontal',ticks=np.linspace(0.00, 50, 6))
    #cb1.set_label('n$_e$ [n$_c$]')
    #cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')

    #cmap = mpl.cm.BrBG
    #Bz = 22.5
    #norm = mpl.colors.Normalize(vmin=-abs(Bz), vmax=abs(Bz))
    #cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap_br,
    #                            norm=norm,
    #                            orientation='horizontal',ticks=np.linspace(-abs(Bz), abs(Bz), 5),alpha=0.7)
    #cb2.set_label(r'E$_y$ [m$_e\omega$/e]')
    #cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')


  ax=plt.subplot(1,2,2)
   
  directory = './txt_1300/'
  px = np.loadtxt(directory+'px2d_x.txt')
  py = np.loadtxt(directory+'py2d_x.txt')
  xx = np.loadtxt(directory+'xx2d_x.txt')
  yy = np.loadtxt(directory+'yy2d_x.txt')
  workx2d = np.loadtxt(directory+'workx2d_x.txt')
  worky2d = np.loadtxt(directory+'worky2d_x.txt')
  fieldex = np.loadtxt(directory+'fieldex2d_x.txt')/4.0
  fieldey = np.loadtxt(directory+'fieldey2d_x.txt')/4.0
  fieldbz = np.loadtxt(directory+'fieldbz2d_x.txt')/4.0
  
  ey_averaged = -8.0/3.2*yy
  bz_averaged = -8.0/3.2*yy
  
  laser_ey = fieldey-ey_averaged
  laser_bz = fieldbz-bz_averaged
  
  
  gg = (px**2+py**2+1)**0.5
  R = gg-px
  theta = np.arctan2(py,px)
  
  number=400
  
  tt = np.linspace(5.0,89.9,850)
  #tt = tt[:,np.newaxis]
  
  select_x = np.array([21,17,2,57,58,68,72,78,172,107])
  
  
      
  line_y1 = np.linspace(-12,-3.2,1001)
  line_x1 = np.zeros_like(line_y1)+5
      
  line_y2 = np.linspace(3.2,12,1001)
  line_x2 = np.zeros_like(line_y2)+5
      
  line_x3 = np.linspace(5,200.0,1001)
  line_y3 = np.zeros_like(line_x3)+3.2
  
  line_x4 = np.linspace(5,200.0,1001)
  line_y4 = np.zeros_like(line_x4)-3.2
  
  for index in select_x:
      norm_x = matplotlib.colors.Normalize(vmin=1.,vmax=501.)
      norm_y = matplotlib.colors.Normalize(vmin=1.,vmax=301.)
  
      plt.scatter(xx[index,50:], yy[index,50:], c=gg[index,50:], norm=norm_x, s=6, cmap='autumn', edgecolors='None')
      plt.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
  
  
  #cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(501), 5) )
  #cbar.set_label('$\gamma$',fontdict=font)
  #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
  #   plt.legend(loc='upper right')
  
  
  
  
  directory = './txt_1300/'
  px = np.loadtxt(directory+'px2d_y.txt')
  py = np.loadtxt(directory+'py2d_y.txt')
  xx = np.loadtxt(directory+'xx2d_y.txt')
  yy = np.loadtxt(directory+'yy2d_y.txt')
  workx2d = np.loadtxt(directory+'workx2d_y.txt')
  worky2d = np.loadtxt(directory+'worky2d_y.txt')
  fieldex = np.loadtxt(directory+'fieldex2d_y.txt')/4.0
  fieldey = np.loadtxt(directory+'fieldey2d_y.txt')/4.0
  fieldbz = np.loadtxt(directory+'fieldbz2d_y.txt')/4.0
  
  ey_averaged = -8.0/3.2*yy
  bz_averaged = -8.0/3.2*yy
  
  laser_ey = fieldey-ey_averaged
  laser_bz = fieldbz-bz_averaged
  
  
  gg = (px**2+py**2+1)**0.5
  R = gg-px
  theta = np.arctan2(py,px)
  
  number=400
  
  tt = np.linspace(5.0,89.9,850)
  #tt = tt[:,np.newaxis]
  
  select_x = np.array([0,1,2,3,12,25,26,29,36,28])
  
  
  line_y1 = np.linspace(-12,-3.2,1001)
  line_x1 = np.zeros_like(line_y1)+5
      
  line_y2 = np.linspace(3.2,12,1001)
  line_x2 = np.zeros_like(line_y2)+5
      
  line_x3 = np.linspace(5,200.0,1001)
  line_y3 = np.zeros_like(line_x3)+3.2
  
  line_x4 = np.linspace(5,200.0,1001)
  line_y4 = np.zeros_like(line_x4)-3.2
  
  for index in select_x:
      norm_x = matplotlib.colors.Normalize(vmin=1.,vmax=501.)
      norm_y = matplotlib.colors.Normalize(vmin=1.,vmax=301.)
  
      plt.scatter(xx[index,50:], yy[index,50:], c=gg[index,50:], norm=norm_y, s=6, cmap='winter', edgecolors='None')
      plt.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
  
  
  #cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(301), 5) )
  #cbar.set_label('$\gamma$',fontdict=font)
  #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
  #   plt.legend(loc='upper right')
  
  plt.xlim(-5,130)
  plt.ylim(-6.0,6.0)
  plt.xlabel('X [$\lambda$]',fontdict=font)
#  plt.ylabel('Y [$\lambda_0$]',fontdict=font)
  plt.xticks(fontsize=25); plt.yticks(fontsize=0);
  #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
  plt.subplots_adjust(left=0.05, bottom=0.15, right=0.98, top=0.97,
                wspace=0.02, hspace=None)


  fig = plt.gcf()
  fig.set_size_inches(24, 6.)
  fig.savefig('./figure_wrap_up/comb_1'+str(n).zfill(4)+'_njp.png',format='png',dpi=160)
  plt.close("all")

