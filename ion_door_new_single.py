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
from scipy.integrate import quad
import scipy.integrate as integrate 
 
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

#  def integrand(x, a, b):
#      return -1/(1.0-1.0/(a*(1-(x/b)**2)+1)**2)**0.5


  charge0 = 6
  mass0   = 1836*12  
  factor  = 0.2 

  for i in [1]:
      fig, ax1 = plt.subplots()
      r0 = 2.7*2*pi
      a0 = np.linspace(1,200,200)

      vy = (charge0*factor*a0*r0/mass0)**0.5
      dt = np.pi/2*(mass0*r0/charge0/factor/a0)**0.5/2/np.pi
      plt.plot(a0,dt,':r',linewidth=4, label='$\delta$t Nonrelativistic',zorder=2)
      dt_array = np.zeros_like(dt)
      for i in range(dt_array.size):
          result_temp = integrate.quad(lambda y: 1.0/(1-1.0/(1+0.5*charge0*factor/r0*a0[i]*(r0**2-y**2)/mass0)**2)**0.5, 0, r0)
          dt_array[i] = result_temp[0]/2/pi
      plt.plot(a0,dt_array,'-r',linewidth=4, label='$\delta$t Relativistic',zorder=2)

      sim_a0=np.array([12.0,27.0,60.0, 120,  190])
      sim_dt=np.array([108.75, 71.75, 52.75, 38.75, 32.45])/3.333*1.05
      #plt.scatter(sim_a0,sim_dt,c='deepskyblue',marker='o',s=200,label='2D-PIC $r_0=3.4\mu m$',edgecolors='black', linewidth=3, alpha=1, zorder=2)
      #sim_a0=np.array([190.0])
      #sim_dt=np.array([10.0])
      plt.scatter(sim_a0,sim_dt,c='tomato',marker='^',s=200,label='$\delta$t 2D-PIC',edgecolors='black', linewidth=3, alpha=1, zorder=3)
      #plt.text(60,19,r'$\Delta t=\sqrt{\frac{A(m_p/m_e)r_0}{2Z\rho a_0}}$'+' w/o Rel',fontdict=font)
      plt.ylim(0,45)
      ax1.tick_params(axis='y', colors='red')
      plt.xlabel('$a_0$',fontdict=font)
      plt.ylabel('$\delta t\ [T_0]$',fontdict=font,color='r')
      plt.xticks([0,50,100,150,200],fontsize=20); 
      plt.yticks([0,10,20,30,40],fontsize=20);
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
      plt.legend(loc='upper center',fontsize=12,framealpha=1.0)
 
      ax2 = ax1.twinx()
      r0 = 2.7*2*pi
      a0 = np.linspace(1,200,200)
      ax2.plot(a0,vy,':b',linewidth=4, label='$v_y$ Nonrelativistic',zorder=1)
      gg = 0.5*charge0*factor*a0*r0/mass0+1.0
      vr = (1-1.0/gg**2)**0.5 
      ax2.plot(a0,vr,'-b',linewidth=4, label='$v_y$ Relativistic',zorder=0)
      sim_a0=np.array([12.0,27.0,60.0, 120,  190])
      sim_vy=np.array([0.114,0.186,0.26,0.33,0.36])
      ax2.scatter(sim_a0,sim_vy,c='deepskyblue',marker='o',s=200,label='$v_y$ 2D-PIC',edgecolors='black', linewidth=3, alpha=1, zorder=2)
      ax2.set_ylabel('$v_y\ [c]$', fontdict=font, color='b')
      ax2.tick_params('y', colors='b')
#      ax2.set_yticklabels(ax2.get_yticklabels(),fontsize=20)
      ax2.yaxis.set_ticks(ticks=[0.0,0.1,0.2,0.3,0.4])
      ax2.yaxis.set_tick_params(labelsize=20)
      ax2.set_ylim(0,0.45)
      #### manifesting colorbar, changing label and axis properties ####
      #plt.ylim(0,0.6)
      #plt.xlabel('$a_0$',fontdict=font)
      #plt.ylabel('$v_y\ [c]$',fontdict=font)
      #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
      plt.legend(loc='lower center',fontsize=12,framealpha=1.0)
      plt.subplots_adjust(left=0.11, bottom=0.13, right=0.87, top=0.95,
                wspace=None, hspace=None)
      #plt.text(40,0.08,r'$v_y=\sqrt{\frac{Z\rho a_0r_0}{A(m_p/m_e)}}$'+' w/o Rel',fontdict=font)



      fig = plt.gcf()
      fig.set_size_inches(7.2, 6.0)
      fig.savefig('./ion_door_new_single.png',format='png',dpi=160)
      plt.close("all")
