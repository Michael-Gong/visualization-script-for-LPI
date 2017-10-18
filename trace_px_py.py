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
#%matplotlib inline
#from colour import Color

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
#from colour import Color
print 'electric field unit: '+str(exunit)
print 'magnetic field unit: '+str(bxunit)
print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 32,  
        }  

######### Parameter you should set ###########
start   =  0  # start time
stop    =  498  # end time
step    =  1  # the interval or step

a0=10
Data = 'Data/'
name = 'electron scattering plot'
total_step = 499
px_2d = np.zeros([3,total_step])
py_2d = np.zeros([3,total_step])
gg_2d = np.zeros([3,total_step])
x_2d  = np.zeros([3,total_step])
######### Script code drawing figure ################

def main(from_path, to_path):
  for n in range(start,stop+step,step):
    #### header data ####
    #plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)
    #ax=plt.subplot(1,1,1,polar=True)
    #fig = plt.figure()
    #ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
    data = sdf.read('Data20/'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    px    = data['Particles/Px/electron'].data/(m0*v0)
    py    = data['Particles/Py/electron'].data/(m0*v0)
    gg    = np.sqrt(px**2+py**2+1.0)
    idi   = data['Particles/ID/electron'].data
    x     = data['Grid/Particles/electron'].data[0]/1.0e-6
    px_2d[0,n-1] = px[np.where(idi == 1)]
    py_2d[0,n-1] = py[np.where(idi == 1)]
    gg_2d[0,n-1] = gg[np.where(idi == 1)]
    x_2d[0,n-1]  = x[np.where(idi == 1)]
    data = sdf.read('Data40/'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    px    = data['Particles/Px/electron'].data/(m0*v0)
    py    = data['Particles/Py/electron'].data/(m0*v0)
    gg    = np.sqrt(px**2+py**2+1.0)
    idi   = data['Particles/ID/electron'].data
    x     = data['Grid/Particles/electron'].data[0]/1.0e-6
    px_2d[1,n-1] = px[np.where(idi == 1)]
    py_2d[1,n-1] = py[np.where(idi == 1)]
    gg_2d[1,n-1] = gg[np.where(idi == 1)]
    x_2d[1,n-1]  = x[np.where(idi == 1)]
    data = sdf.read('Data80/'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    px    = data['Particles/Px/electron'].data/(m0*v0)
    py    = data['Particles/Py/electron'].data/(m0*v0)
    gg    = np.sqrt(px**2+py**2+1.0)
    idi   = data['Particles/ID/electron'].data
    x     = data['Grid/Particles/electron'].data[0]/1.0e-6
    px_2d[2,n-1] = px[np.where(idi == 1)]
    py_2d[2,n-1] = py[np.where(idi == 1)]
    gg_2d[2,n-1] = gg[np.where(idi == 1)]
    x_2d[2,n-1]  = x[np.where(idi == 1)]
    print 'finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%'
  
  plt.subplot(2,1,1)
  x = np.linspace(-1.0,1.0,501)
  y = x**2
  plt.plot(x,y,'-k',linewidth=3,label='analytical solution')
  a0=20
  plt.scatter(py_2d[0,:]/a0,px_2d[0,:]/(a0**2)*2,s=0.5,c=(0.0,0.0,192.0/255.0),label=r'dx=1/20$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.scatter(py_2d[0,:]/a0,px_2d[0,:]/(a0**2)*2,s=0.5,c=(0.0,192.0/255.0,0.0),label=r'dx=1/40$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.scatter(py_2d[0,:]/a0,px_2d[0,:]/(a0**2)*2,s=0.5,c=(192.0/255.0,0.0,0.0),label=r'dx=1/80$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.xlabel('$p_y/a_0$',fontdict=font)
  plt.ylabel('$2p_x/a_0^2$',fontdict=font)
  plt.xticks(fontsize=18) 
  plt.yticks(fontsize=18)
  plt.xlim(-1.5,1.5)
  plt.ylim(-0.2,1.2)
  #ax.set_rticks(20);
  plt.title(r'$p_x-p_y$',fontdict=font)
  #ax.set_title(r"(a) $\xi_0=250$", va='bottom', fontdict=font)
  plt.legend(loc='upper center',framealpha=0.0,markerscale=10.0,fontsize=20.0,scatterpoints=3)

  plt.subplot(2,1,2)
  x = np.linspace(0.0,100.0,501)
  y = x*0.0+1
  plt.plot(x,y,'-k',linewidth=3,label='analytical solution')
  plt.scatter(x_2d[0,:],gg_2d[0,:]-px_2d[0,:],s=2.0,c=(0.0,0.0,192.0/255.0),label=r'dx=1/20$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.scatter(x_2d[0,:],gg_2d[0,:]-px_2d[0,:],s=2.0,c=(0.0,192.0/255.0,0.0),label=r'dx=1/40$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.scatter(x_2d[0,:],gg_2d[0,:]-px_2d[0,:],s=2.0,c=(192.0/255.0,0.0,0.0),label=r'dx=1/80$\lambda; a_0=$'+str(a0),edgecolors='None')
  plt.xlabel('$x [\lambda]$',fontdict=font)
  plt.ylabel('$\gamma-p_x$',fontdict=font)
  plt.xticks(fontsize=18) 
  plt.yticks(fontsize=18)
  plt.xlim(0,100)
  plt.ylim(0.8,4.0)
  #ax.set_rticks(20);
  plt.title(r'$\gamma-px\ vs\ x$',fontdict=font)
  #ax.set_title(r"(a) $\xi_0=250$", va='bottom', fontdict=font)
  #plt.legend(loc='upper center',framealpha=0.0,markerscale=10.0,fontsize=20.0,scatterpoints=3)
  fig = plt.gcf()
  fig.set_size_inches(10,16)
  fig.savefig('px_py'+'.png',format='png',dpi=100)
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
         #       exit()
        print "from path:", option.from_path
        print "to path:", option.to_path
        if not os.path.exists(option.to_path):
                os.mkdir(option.to_path)
        main(option.from_path,option.to_path)

 
