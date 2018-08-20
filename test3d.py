
import sdf 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal



######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  



######### Parameter you should set ###########
start   =  1  # start time
stop    =  88  # end time
step    =  1  # the interval or step

#youwant = ['']
youwant =  ['ey','ex','Subset_high_e_density','Subset_high_e_ekbar','ey_averaged','bz_averaged']
#youwant.append('carbon_ekbar')
#youwant.append('positron_ekbar')
#youwant.append('electron_en')
#youwant.append('photon_en')
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...
#youwant dist_fn electron_x_px...

from_path = './Data_3d/'
to_path   = './Data_3d/'


######### Script code drawing figure ################
for n in range(start,stop+step,step):
  #### header data ####
  data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  X, Y = np.meshgrid(x, y)
  
  for name in youwant:
    if (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
              eexx = data['Electric Field/'+str.capitalize(name)].data/exunit
              n3d=len(eexx[0,0,:])
              ex = (eexx[:,:,n3d//2-1]+eexx[:,:,n3d//2])/2 
              if np.min(ex.T) == np.max(ex.T):
                  continue
              eee=np.max([-np.min(ex.T),np.max(ex.T)])
              levels = np.linspace(-eee, eee, 40)
              plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdBu)
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
              cbar.set_label('Normalized electric field',fontdict=font)        
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Y [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              fig = plt.gcf()
              fig.set_size_inches(12, 7)
              fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
    elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
              eexx = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
              n3d=len(eexx[0,0,:])
              ex = (eexx[:,:,n3d//2-1]+eexx[:,:,n3d//2])/2
              if np.min(ex.T) == np.max(ex.T):
                  continue
              eee=np.max([-np.min(ex.T),np.max(ex.T)])
              levels = np.linspace(-eee, eee, 40)
              plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdBu)
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
              cbar.set_label('Normalized magnetic field',fontdict=font)        
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Y [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              fig = plt.gcf()
              fig.set_size_inches(12, 7)
              fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
    elif (name[-7:] == 'density'):
              ddeen = data['Derived/Number_Density/'+name[0:-8]].data/denunit
              n3d=len(ddeen[0,0,:])
              den = (ddeen[:,:,n3d//2-1]+ddeen[:,:,n3d//2])/2
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
              ddeen = data['Derived/EkBar/'+name[0:-6]].data/(q0*1.0e6)
              n3d=len(ddeen[0,0,:])
              den = (ddeen[:,:,n3d//2-1]+ddeen[:,:,n3d//2])/2
              if np.min(den.T) == np.max(den.T):
                  continue
              levels = np.linspace(np.min(den.T), np.max(den.T), 40) 
              plt.contourf(X, Y, den.T, levels=levels, cmap=cm.nipy_spectral)
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
    elif (name[-8:] == 'theta_en'):
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
              den = data['dist_fn/en/'+name[0:-3]].data[:,0,0]
              dist_x  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)
              plt.plot(dist_x,den,'-r',linewidth=3)
              #### manifesting colorbar, changing label and axis properties ####
              plt.xlabel('Energy [MeV]',fontdict=font)
              plt.ylabel('dN/dE [A.U.]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              fig = plt.gcf()
              fig.set_size_inches(12, 7)
              fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
  print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

