#!/usr/local/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os

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
print 'electric field unit: '+str(exunit)
print 'magnetic field unit: '+str(bxunit)
print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  



######### Parameter you should set ###########
start   =  1  # start time
stop    =  59  # end time
step    =  1  # the interval or step
youwant = ['ey','electron_density','electron_ekbar']
#youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey']
#youwant =  ['bz','ex','ey_averaged','ez','electron_density','carbon_density','photon_density','positron_density','electron_ekbar','photon_ekbar','electron_x_px']
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...
#youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...

######### Script code drawing figure ################
def main(from_path,to_path):
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    for name in youwant:
      if (name[0:2] == 'ex') or (name[0:2] == 'ey') or (name[0:2] == 'ez'):
                ex = data['Electric Field/'+str.capitalize(name)].data/exunit
                if np.min(ex.T) == np.max(ex.T):
                    continue
                plt.plot(x,ex,'-r',linewidth=3)
                #### manifesting colorbar, changing label and axis properties ####
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Normalized electric field',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
      elif (name[0:2] == 'bx') or (name[0:2] == 'by') or (name[0:2] == 'bz'):
                ex = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
                if np.min(ex.T) == np.max(ex.T):
                    continue
                plt.plot(x,ex,'-r',linewidth=3)
                #### manifesting colorbar, changing label and axis properties ####
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel('Normalized magnetic field',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
      elif (name[-7:] == 'density'):
                den = data['Derived/Number_Density/'+name[0:-8]].data/denunit
                if np.min(den.T) == np.max(den.T):
                    continue
                plt.plot(x,den,'-r',linewidth=3)
                #### manifesting colorbar, changing label and axis properties ####
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel(name+'[$n_c$]',fontdict=font)
                plt.xticks(fontsize=20); plt.yticks(fontsize=20);
                plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
                fig = plt.gcf()
                fig.set_size_inches(12, 7)
                fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
                plt.close("all")
      elif (name[-5:] == 'ekbar'):
                den = data['Derived/EkBar/'+name[0:-6]].data/(q0*1.0e6)
                if np.min(den.T) == np.max(den.T):
                    continue
                plt.plot(x,den,'-r',linewidth=3)
                #### manifesting colorbar, changing label and axis properties ####
                plt.xlabel('X [$\mu m$]',fontdict=font)
                plt.ylabel(name+'[MeV]',fontdict=font)
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
                den = data['dist_fn/theta_en/'+name[0:-9]].data[:,:,0]
                den = np.log(den+1.0)
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
    print 'finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%'
  
  
if __name__ == "__main__":
  	parser = OptionParser()
  	parser.add_option("-f","--from_path",
  		dest = "from_path",
  		type = "string",
  		default = "Data")
  	parser.add_option("-t","--to_path",
  		dest = "to_path",
  		type = "string",
  		default = "jpg")
  	(option,args) = parser.parse_args()
  	if option.from_path[-1:] != '/' :
  		option.from_path += '/'
  	option.to_path = option.to_path 
  	if option.to_path[-1:] != '/' :
  		option.to_path += '/'
  	if not os.path.exists(option.from_path):
  		print 'error: input data path not exist'
  		exit()
  	print "from path:", option.from_path
  	print "to path:", option.to_path
	if not os.path.exists(option.to_path):
		os.mkdir(option.to_path)
	main(option.from_path, option.to_path)
