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

import multiprocessing as mp
 


 
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



def processplot(n):
  from_path_list=['./a0_200_n_02_l_10_qe/','./a0_200_n_02_l_10_no/','./a0_100_n_01_l_10_qe/','./a0_100_n_01_l_10_no/']
  to_path_list=['./jpg_a0_200_n_02_l_10_qe/','./jpg_a0_200_n_02_l_10_no/','./jpg_a0_100_n_01_l_10_qe/','./jpg_a0_100_n_01_l_10_no/'] 
  
  ######### Parameter you should set ###########
  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  youwant =  ['jx_averaged','jy_averaged','ey','bz','ex','ey_averaged','bz_averaged','e_0_density','e_0_ekbar','e_1_density','e_1_ekbar','c_0_density','c_0_ekbar','e_0_en','e_1_en','c_0_en']#,'e_0_theta_en','e_1_theta_en','c_0_theta_en','c_1_theta_en','e_0_x_px','e_1_x_px','c_0_x_px','c_1_x_px']
#  youwant = ['e_0_en','e_1_en','c_0_en']#,'e_0_theta_en','e_1_theta_en','c_0_theta_en','c_1_theta_en','e_0_x_px','e_1_x_px','c_0_x_px','c_1_x_px']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  ######### Script code drawing figure ################
  for name in youwant:
    for i in range(np.size(from_path_list)):
        from_path = from_path_list[i]
        to_path   = to_path_list[i]
        if (name[0:2] == 'jx') or (name[0:2] == 'jy') or (name[0:2] == 'jz'):
                data = sdf.read(from_path+"current"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                x  = data['Grid/Grid_mid'].data[0]/1.0e-6
                print('ok')
                y  = data['Grid/Grid_mid'].data[1]/1.0e-6
                X, Y = np.meshgrid(x, y)
                ex = data['Current/'+str.capitalize(name)].data/jalf
                if np.min(ex.T) == np.max(ex.T):
                    continue
                eee=np.max([-np.min(ex.T),np.max(ex.T)])
                levels = np.linspace(-eee, eee, 40)
                plt.contourf(X, Y, ex.T, norm=MidpointNormalize(midpoint=0.), cmap=cm.PiYG)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar()
                cbar.set_label('Normalized current',fontdict=font)
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Y [$\mu m$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                #plt1 = plt.twinx()
                #plt1.plot(x,ex[:,y.size/2.0],'-k',linewidth=2.5)
                #plt1.set_ylabel('Normalized '+name)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
                
        if (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
                data = sdf.read(from_path+"e_fields"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                x  = data['Grid/Grid_mid'].data[0]/1.0e-6
                print('ok')
                y  = data['Grid/Grid_mid'].data[1]/1.0e-6
                X, Y = np.meshgrid(x, y)
                ex = data['Electric Field/'+str.capitalize(name)].data/exunit
                if np.min(ex.T) == np.max(ex.T):
                    continue
                eee=np.max([-np.min(ex.T),np.max(ex.T)])
                levels = np.linspace(-eee, eee, 40)
                plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.seismic)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
                cbar.set_label('Normalized electric field',fontdict=font)        
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Y [$\mu m$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                #plt1 = plt.twinx()
                #plt1.plot(x,ex[:,y.size/2.0],'-k',linewidth=2.5)
                #plt1.set_ylabel('Normalized '+name)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")

        elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
                data = sdf.read(from_path+"b_fields"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                x  = data['Grid/Grid_mid'].data[0]/1.0e-6
                print('ok')
                y  = data['Grid/Grid_mid'].data[1]/1.0e-6
                X, Y = np.meshgrid(x, y)
                ex = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
                if np.min(ex.T) == np.max(ex.T):
                    continue
                eee=np.max([-np.min(ex.T),np.max(ex.T)])
                levels = np.linspace(-eee, eee, 40)
                plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.seismic)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
                cbar.set_label('Normalized magnetic field',fontdict=font)        
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Y [$\mu m$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                #plt1 = plt.twinx()
                #plt1.plot(x,ex[:,y.size/2.0],'-k',linewidth=2.5)
                #plt1.set_ylabel('Normalized '+name)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")

        elif (name[-7:] == 'density'):
                data = sdf.read(from_path+"density"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                x  = data['Grid/Grid_mid'].data[0]/1.0e-6
                print('ok')
                y  = data['Grid/Grid_mid'].data[1]/1.0e-6
                X, Y = np.meshgrid(x, y)
                den = data['Derived/Number_Density/'+name[0:-8]].data/denunit
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40) 
                plt.contourf(X, Y, den.T, levels=levels, cmap=cm.nipy_spectral)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[$n_c$]', fontdict=font)
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Y [$\mu m$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")

        elif (name[-5:] == 'ekbar'):
                data = sdf.read(from_path+"ekbar"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                x  = data['Grid/Grid_mid'].data[0]/1.0e-6
                print('ok')
                y  = data['Grid/Grid_mid'].data[1]/1.0e-6
                X, Y = np.meshgrid(x, y)
                den = data['Derived/EkBar_averaged/'+name[0:-6]].data/(q0*1.0e6)
                if np.min(den.T) == np.max(den.T):
                    continue
                levels = np.linspace(np.min(den.T), np.max(den.T), 40) 
                plt.contourf(X, Y, den.T, levels=levels, cmap=cm.jet)
                #### manifesting colorbar, changing label and axis properties ####
                cbar=plt.colorbar(ticks=np.linspace(np.min(den.T), np.max(den.T), 5))
                cbar.set_label(name+'[MeV]', fontdict=font)
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Y [$\mu m$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
        elif (name[-4:] == 'x_px'):
                data = sdf.read(from_path+"dist"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
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
                data = sdf.read(from_path+"dist"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
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
                data = sdf.read(from_path+"dist"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
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
                data = sdf.read(from_path+"dist"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                denden = data['dist_fn/theta_en/'+name[0:-9]].data[:,:,0]
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
                data = sdf.read(from_path+"dist"+str(n).zfill(4)+".sdf",dict=True)
                header=data['Header']
                time=header['time']
                den = data['dist_fn/en/'+name[0:-3]].data[:,0,0]
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

if __name__ == '__main__':
  start   = 0  # start time
  stop    = 12 # end time
  step    = 1

  inputs = range(start,stop+step,step)
  pool   = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)
