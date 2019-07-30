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
  print('current density unit nc: '+str(jalf))
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 18,  
          }  
  font_size   = 20
  font_size2  = 18
  color_level = 49
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####

  def counter_derived_variable(x,y,value,name,c_max,c_min,c_map,scale,time):
      if c_max == c_min:
	      return 0
      if scale == 'log':
          levels = np.logspace(np.log10(c_min),np.log10(c_max),color_level)
          plt.contourf(X, Y, value.T, levels=levels, norm=colors.LogNorm(vmin=c_min, vmax=c_max))
          cbar=plt.colorbar(pad=0.01, ticks=np.logspace(np.log10(c_min), np.log10(c_max), 5))
      if scale == 'linear':
          levels = np.linspace(c_min,c_max,color_level)
          plt.contourf(X, Y, value.T, levels=levels, norm=colors.Normalize(vmin=c_min, vmax=c_max))
          cbar=plt.colorbar(pad=0.01, ticks=np.linspace(c_min, c_max, 5))
      if scale == 'bilinear':
          levels = np.linspace(c_min,c_max,color_level)    
          plt.contourf(X, Y, value.T, levels=levels, norm=MidpointNormalize(midpoint=0.), cmap=c_map)
          cbar=plt.colorbar(pad=0.01, ticks=np.linspace(c_min, c_max, 5))
      #### manifesting colorbar, changing label and axis properties ####
      if (name[0:2] == 'jx') or (name[0:2] == 'jy') or (name[0:2] == 'jz'):
          cbar.set_label('Normalized current density '+r'($\alpha=j_0\lambda^2/4\pi J_A$)',fontdict=font2)
      elif (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
          cbar.set_label('Normalized electric field',fontdict=font2)
      elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
          cbar.set_label('Normalized magnetic field',fontdict=font2)
      elif (name[-10:-1] == 'polarize_'):
          cbar.set_label('Polarization',fontdict=font2)
      elif (name[-7:] == 'density'):
          cbar.set_label(name+' $[n_c]$',fontdict=font2)
      elif (name[-5:] == 'ekbar'):
          cbar.set_label(name+' [MeV]',fontdict=font2)
      cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])
      plt.xlabel('x [$\mu m$]',fontdict=font); plt.ylabel('y [$\mu m$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(9.8, 7)
      fig.savefig('./Data/'+name+str(n).zfill(4)+'.png',format='png',dpi=80)
      plt.close("all")

 ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  14  # end time
  step    =  1  # the interval or step
  from_path='./Data/' 
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  #youwant =  ['electron_en','electron_no_en','ey','ex','ey_averaged','bz','bz_averaged','electron_density','electron_ekbar']#,'Subset_high_e_density','Subset_high_e_ekbar']
  youwant =  ['ex','ey','ey_averaged','bz','bz_averaged','electron_s_density','electron_s_ekbar','ion_s_density','ion_s_ekbar','electron_s_polarize_x','electron_s_polarize_y','electron_s_polarize_z','ion_s_polarize_x','ion_s_polarize_y','ion_s_polarize_z']
#  youwant =  ['ex','ey','ey_averaged','bz','bz_averaged','Electron_density','Electron_ekbar','photon_density','photon_ekbar','Electron_theta_en','Electron_en','Carbon_theta_en','Carbon_en','photon_theta_en','photon_en']#,'Subset_high_e_density','Subset_high_e_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y)
    
    for name in youwant:
        if (name[0:2] == 'jx') or (name[0:2] == 'jy') or (name[0:2] == 'jz'):
            jx = data['Current/'+str.capitalize(name)].data/jalf
            mmm = max(np.max(jx),-np.min(jx))
            counter_derived_variable(x=X,y=Y,value=jx, name=name, c_max=mmm, c_min=-mmm, c_map='PiYG', scale='bilinear', time=time)
        elif (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
            ex = data['Electric Field/'+str.capitalize(name)].data/exunit
            mmm = max(np.max(ex),-np.min(ex))
            counter_derived_variable(x=X,y=Y,value=ex, name=name, c_max=mmm, c_min=-mmm, c_map='RdBu_r', scale='bilinear', time=time)
        elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
            bx = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
            mmm = max(np.max(bx),-np.min(bx))
            counter_derived_variable(x=X,y=Y,value=bx, name=name, c_max=mmm, c_min=-mmm, c_map='BrBG_r', scale='bilinear', time=time)
        elif (name[-10:-1] == 'polarize_'):
            pol = data['Derived/Average_Particle_spin_'+name[-1]+'/'+name[:-11]].data
            mmm = max(np.max(pol),-np.min(pol))
            counter_derived_variable(x=X,y=Y,value=pol, name=name, c_max=mmm, c_min=-mmm, c_map='bwr', scale='bilinear', time=time)
        elif (name[-7:] == 'density'):
            den = data['Derived/Number_Density/'+name[0:-8]].data/denunit
            counter_derived_variable(x=X,y=Y,value=den, name=name, c_max=np.max(den), c_min=np.min(den), c_map='viridis', scale='linear', time=time)
        elif (name[-5:] == 'ekbar'):
            ek = data['Derived/Average_Particle_Energy/'+name[0:-6]].data/(q0*1.0e6)
            counter_derived_variable(x=X,y=Y,value=ek, name=name, c_max=np.max(ek), c_min=np.min(ek), c_map='magma', scale='linear', time=time)
        elif (name[-4:] == 'x_px'):
                den = data['dist_fn/x_px/'+name[0:-5]].data[:,:,0]
                den = np.log(den+1.0)
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40)
                dist_x  = data['Grid/x_px/'+name[0:-5]].data[0]/1.0e-6
                dist_y  = data['Grid/x_px/'+name[0:-5]].data[1]/(m0*v0)
                dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
                plt.contourf(dist_X, dist_Y, den.T, levels=levels, cmap=cm.nipy_spectral)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[$log_{10}(A.U.)$]', fontdict=font)
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('$P_x$ [$m_ec$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
        elif (name[-4:] == 'y_py'):
                den = data['dist_fn/y_py/'+name[0:-5]].data[:,:,0]
                den = np.log(den+1.0)
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40)
                dist_x  = data['Grid/y_py/'+name[0:-5]].data[0]/1.0e-6
                dist_y  = data['Grid/y_py/'+name[0:-5]].data[1]/(m0*v0)
                dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
                plt.contourf(dist_X, dist_Y, den.T, levels=levels, cmap=cm.nipy_spectral)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[$log_{10}(A.U.)$]', fontdict=font)
                plt.xlabel('Y [$\mu m$]',fontdict=font)
                plt.ylabel('$P_y$ [$m_ec$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
        elif (name[-5:] == 'py_pz'):
                den = data['dist_fn/py_pz/'+name[0:-6]].data[:,:,0]
                den = np.log(den+1.0)
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40)
                dist_x  = data['Grid/py_pz/'+name[0:-6]].data[0]/(m0*v0)
                dist_y  = data['Grid/py_pz/'+name[0:-6]].data[1]/(m0*v0)
                dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
                plt.contourf(dist_X, dist_Y, den.T, levels=levels, cmap=cm.nipy_spectral)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[$log_{10}(A.U.)$]', fontdict=font)
                plt.xlabel('$P_y$ [$m_ec$]',fontdict=font)
                plt.ylabel('$P_z$ [$m_ec$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
        elif (name[-8:] == 'theta_en'):
                denden = data['dist_fn/theta_en/'+name[0:-9]].data
                den = np.log(denden+1.0)
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40)
                dist_x  = data['Grid/theta_en/'+name[0:-9]].data[0]
                dist_y  = data['Grid/theta_en/'+name[0:-9]].data[1]/(q0*1.0e6)
                dist_X, dist_Y = np.meshgrid(dist_x, dist_y)
                plt.contourf(dist_X, dist_Y, den.T, levels=levels, cmap=cm.nipy_spectral)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[$log_{10}(A.U.)$]', fontdict=font)
                plt.xlabel('$\Psi$ [rad]',fontdict=font)
                plt.ylabel('$Energy$ [MeV]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                plt1 = plt.twinx()
                plt1.plot(dist_x,np.sum(denden,axis=1),'-y',linewidth=2.5)
                #plt1.set_ylabel('Normalized '+name)  
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
        elif (name[-2:] == 'en'):
                den = data['dist_fn/en/'+name[0:-3]].data
                dist_x  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)
                plt.plot(dist_x,den,'-r',linewidth=3)
                #### manifesting colorbar, changing label and axis properties ####
                plt.xlabel('Energy [MeV]',fontdict=font)
                plt.ylabel('dN/dE [A.U.]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.yscale('log')
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  
