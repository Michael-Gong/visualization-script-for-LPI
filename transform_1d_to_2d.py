from scipy.integrate import odeint
#%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

font = {'family' : 'helvetica',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
        }

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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

######### Parameter you should set ###########
start   =  1  # start time
stop    =  399  # end time
step    =  1  # the interval or step


youwant = ['ey','electron_density','electron_ekbar','bz']
#youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey']
#youwant =  ['bz','ex','ey_averaged','ez','electron_density','carbon_density','photon_density','positron_density','electron_ekbar','photon_ekbar','electron_x_px']
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...

nx=10000
nt=400

x_grid=np.array(nx)
t_grid=np.array(nt)

ey      = np.zeros([nx,nt])
eden    = np.zeros([nx,nt])
ekbar   = np.zeros([nx,nt])
bz      = np.zeros([nx,nt]) 


######### Script code drawing figure ################
def main(from_path,to_path):
  for n in range(start,stop,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    t_grid[n-start]=header['time']*v0/wavelength
    x_grid=data['Grid/Grid_mid'].data[0]/wavelength
    for name in youwant:
      if (name == 'ey'):
                ey[:,n-start] = data['Electric Field/'+str.capitalize(name)].data/exunit
      elif (name == 'bz'):
                bz[:,n-start] = data['Magnetic Field/'+str.capitalize(name)].data/bxunit
      elif (name == 'electron_density'):
                eden[:,n-start] = data['Derived/Number_Density/'+name[0:-8]].data/denunit
      elif (name == 'electron_ekbar'):
                ekbar[:,n-start] = data['Derived/EkBar/'+name[0:-6]].data/(q0*1.0e6)
    print 'finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%'




  ax=plt.subplot(1,1,1)
  #### manifesting colorbar, changing label and axis properties ####
  eee=np.max(np.abs(ey))
  levels = np.linspace(-eee, eee, 32)   
  plt.contourf(x_grid, t_grid, ey, levels=levels, cmap=cm.hot, alpha=1.0)
  cbar=plt.colorbar(ticks=np.linspace(-eee, eee, 5))
  cbar.set_label('Normalized electric field',fontdict=font)
  #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
  plt.xlabel(r'$x\ [\lambda]$',fontdict=font)
  plt.ylabel(r'$t\ [T_0]$',fontdict=font)
  plt.xticks(fontsize=26.0); plt.yticks(fontsize=26);
  #plt.ylim(-1.025,1.025)
  #plt.xlim(start_x,start_x+length_x)
  #plt.legend(loc='best')
  #plt.show()
  plt.subplots_adjust(top=0.92, bottom=0.08, left=0.15, right=0.95, hspace=0.05, wspace=0.30)
  fig = plt.gcf()
  fig.set_size_inches(8, 6.5)
  #fig.set_size_inches(5, 4.5)
  fig.savefig('./ey_1d_2d.png',format='png',dpi=160)
  plt.close("all")




  ax=plt.subplot(1,1,1)
  #### manifesting colorbar, changing label and axis properties ####
  levels = np.linspace(np.min(eden), np.max(eden), 32)   
  plt.contourf(x_grid, t_grid, eden, levels=levels, cmap=cm.bwr, alpha=1.0)
  cbar=plt.colorbar(ticks=np.linspace(np.min(eden), np.max(eden), 5))
  cbar.set_label('Electron density',fontdict=font)
  #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
  plt.xlabel(r'$x\ [\lambda]$',fontdict=font)
  plt.ylabel(r'$t\ [T_0]$',fontdict=font)
  plt.xticks(fontsize=26.0); plt.yticks(fontsize=26);
  #plt.ylim(-1.025,1.025)
  #plt.xlim(start_x,start_x+length_x)
  #plt.legend(loc='best')
  #plt.show()
  plt.subplots_adjust(top=0.92, bottom=0.08, left=0.15, right=0.95, hspace=0.05, wspace=0.30)
  fig = plt.gcf()
  fig.set_size_inches(8, 6.5)
  #fig.set_size_inches(5, 4.5)
  fig.savefig('./eden_1d_2d.png',format='png',dpi=160)
  plt.close("all")




  ax=plt.subplot(1,1,1)
  #### manifesting colorbar, changing label and axis properties ####
  levels = np.linspace(np.min(ekbar), np.max(ekbar), 32)   
  plt.contourf(x_grid, t_grid, ekbar, levels=levels, cmap=cm.magma, alpha=1.0)
  cbar=plt.colorbar(ticks=np.linspace(np.min(ekbar), np.max(ekbar), 5))
  cbar.set_label('Electron ekbar',fontdict=font)
  #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
  plt.xlabel(r'$x\ [\lambda]$',fontdict=font)
  plt.ylabel(r'$t\ [T_0]$',fontdict=font)
  plt.xticks(fontsize=26.0); plt.yticks(fontsize=26);
  #plt.ylim(-1.025,1.025)
  #plt.xlim(start_x,start_x+length_x)
  #plt.legend(loc='best')
  #plt.show()
  plt.subplots_adjust(top=0.92, bottom=0.08, left=0.15, right=0.95, hspace=0.05, wspace=0.30)
  fig = plt.gcf()
  fig.set_size_inches(8, 6.5)
  #fig.set_size_inches(5, 4.5)
  fig.savefig('./ekbar_1d_2d.png',format='png',dpi=160)
  plt.close("all")




  ax=plt.subplot(1,1,1)
  #### manifesting colorbar, changing label and axis properties ####
  eee=np.max(np.abs(bz))
  levels = np.linspace(-eee, eee, 32)   
  plt.contourf(x_grid, t_grid, bz, levels=levels, cmap=cm.hot, alpha=1.0)
  cbar=plt.colorbar(ticks=np.linspace(-eee, eee, 5))
  cbar.set_label('Normalized magnetic field',fontdict=font)
  #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
  plt.xlabel(r'$x\ [\lambda]$',fontdict=font)
  plt.ylabel(r'$t\ [T_0]$',fontdict=font)
  plt.xticks(fontsize=26.0); plt.yticks(fontsize=26);
  #plt.ylim(-1.025,1.025)
  #plt.xlim(start_x,start_x+length_x)
  #plt.legend(loc='best')
  #plt.show()
  plt.subplots_adjust(top=0.92, bottom=0.08, left=0.15, right=0.95, hspace=0.05, wspace=0.30)
  fig = plt.gcf()
  fig.set_size_inches(8, 6.5)
  #fig.set_size_inches(5, 4.5)
  fig.savefig('./bz_1d_2d.png',format='png',dpi=160)
  plt.close("all")




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
