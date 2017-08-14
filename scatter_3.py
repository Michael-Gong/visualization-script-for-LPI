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
mport matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
#from colour import Color
print 'electric field unit: '+str(exunit)
print 'magnetic field unit: '+str(bxunit)
print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  



######### Parameter you should set ###########
start   =  0  # start time
stop    =  211  # end time
step    =  1  # the interval or step

Data = 'Data0/'
name = 'electron scattering plot'
######### Script code drawing figure ################

def main(from_path, to_path):
  for n in range(start,stop+step,step):
    #### header data ####
    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,wspace=0.25,hspace=0.3)

    plt.subplot(1,3,1)
    data0 = sdf.read("./epoch2dno/Data0/"+str(n).zfill(4)+".sdf",dict=True)
    data1 = sdf.read("./epoch2drr/Data0/"+str(n).zfill(4)+".sdf",dict=True)
    data2 = sdf.read("./epoch2dqe/Data0/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    x1  = data1['Grid/Particles/electron'].data[0]/1.0e-6
    y1  = data1['Grid/Particles/electron'].data[1]/1.0e-6
    x2  = data2['Grid/Particles/electron'].data[0]/1.0e-6
    y2  = data2['Grid/Particles/electron'].data[1]/1.0e-6
    x   = data0['Grid/Grid_mid'].data[0]/1.0e-6
    y   = data0['Grid/Grid_mid'].data[1]/1.0e-6
    X,Y = np.meshgrid(x,y)
    ex  = data0['Electric Field/Ey'].data/exunit
    if np.min(ex.T) == np.max(ex.T):
          continue
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-eee, eee, 24)
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdGy)
  #  plt.scatter(x0,y0,s=10,c=Color(rgb=(1,0,0)),label='1',edgecolors='None') 
  #  plt.scatter(x1,y1,s=10,c=Color(rgb=(0,1,0)),label='1',edgecolors='None') 
  #  plt.scatter(x2,y2,s=10,c=Color(rgb=(0,0,1)),label='1',edgecolors='None') 
    plt.scatter(x2,y2,s=8,c=(192.0/255.0,0,0),label='QED RR',edgecolors='None') 
    plt.scatter(x1,y1,s=8,c=(0,192.0/255.0,0),label='LL RR',edgecolors='None') 
    plt.scatter(x0,y0,s=8,c=(0,0,225.0/255.0),label='No RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=4.0,fontsize=20.0)
    plt.xlim(0,100)
    plt.ylim(-50,50)
    plt.text(5,45,r'$\xi_0=250$',fontsize=20)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Y [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.title(name+' at '+str(round(time/3.3333e-15,6))+' $T_0$',fontdict=font)

    plt.subplot(1,3,2)
    data0 = sdf.read("./epoch2dno/Data1/"+str(n).zfill(4)+".sdf",dict=True)
    data1 = sdf.read("./epoch2drr/Data1/"+str(n).zfill(4)+".sdf",dict=True)
    data2 = sdf.read("./epoch2dqe/Data1/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    x1  = data1['Grid/Particles/electron'].data[0]/1.0e-6
    y1  = data1['Grid/Particles/electron'].data[1]/1.0e-6
    x2  = data2['Grid/Particles/electron'].data[0]/1.0e-6
    y2  = data2['Grid/Particles/electron'].data[1]/1.0e-6
    x   = data0['Grid/Grid_mid'].data[0]/1.0e-6
    y   = data0['Grid/Grid_mid'].data[1]/1.0e-6
    X,Y = np.meshgrid(x,y)
    ex  = data0['Electric Field/Ey'].data/exunit
    if np.min(ex.T) == np.max(ex.T):
          continue
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-eee, eee, 24)
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdGy)
  #  plt.scatter(x0,y0,s=10,c=Color(rgb=(1,0,0)),label='1',edgecolors='None') 
  #  plt.scatter(x1,y1,s=10,c=Color(rgb=(0,1,0)),label='1',edgecolors='None') 
  #  plt.scatter(x2,y2,s=10,c=Color(rgb=(0,0,1)),label='1',edgecolors='None') 
    plt.scatter(x2,y2,s=8,c=(192.0/255.0,0,0),label='QED RR',edgecolors='None') 
    plt.scatter(x1,y1,s=8,c=(0,192.0/255.0,0),label='LL RR',edgecolors='None') 
    plt.scatter(x0,y0,s=8,c=(0,0,225.0/255.0),label='No RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=4.0,fontsize=20.0)
    plt.xlim(0,100)
    plt.ylim(-50,50)
    plt.text(5,45,r'$\xi_0=350$',fontsize=20)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Y [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.title(name+' at '+str(round(time/3.3333e-15,6))+' $T_0$',fontdict=font)

    plt.subplot(1,3,3)
    data0 = sdf.read("./epoch2dno/Data2/"+str(n).zfill(4)+".sdf",dict=True)
    data1 = sdf.read("./epoch2drr/Data2/"+str(n).zfill(4)+".sdf",dict=True)
    data2 = sdf.read("./epoch2dqe/Data2/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    x1  = data1['Grid/Particles/electron'].data[0]/1.0e-6
    y1  = data1['Grid/Particles/electron'].data[1]/1.0e-6
    x2  = data2['Grid/Particles/electron'].data[0]/1.0e-6
    y2  = data2['Grid/Particles/electron'].data[1]/1.0e-6
    x   = data0['Grid/Grid_mid'].data[0]/1.0e-6
    y   = data0['Grid/Grid_mid'].data[1]/1.0e-6
    X,Y = np.meshgrid(x,y)
    ex  = data0['Electric Field/Ey'].data/exunit
    if np.min(ex.T) == np.max(ex.T):
          continue
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-eee, eee, 24)
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdGy)
  #  plt.scatter(x0,y0,s=10,c=Color(rgb=(1,0,0)),label='1',edgecolors='None') 
  #  plt.scatter(x1,y1,s=10,c=Color(rgb=(0,1,0)),label='1',edgecolors='None') 
  #  plt.scatter(x2,y2,s=10,c=Color(rgb=(0,0,1)),label='1',edgecolors='None') 
    plt.scatter(x2,y2,s=8,c=(192.0/255.0,0,0),label='QED RR',edgecolors='None') 
    plt.scatter(x1,y1,s=8,c=(0,192.0/255.0,0),label='LL RR',edgecolors='None') 
    plt.scatter(x0,y0,s=8,c=(0,0,225.0/255.0),label='No RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=4.0,fontsize=20.0)
    plt.xlim(0,100)
    plt.ylim(-50,50)
    plt.text(5,45,r'$\xi_0=450$',fontsize=20)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Y [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.title(name+' at '+str(round(time/3.3333e-15,6))+' $T_0$',fontdict=font)

    fig = plt.gcf()
    fig.set_size_inches(30, 9)
    fig.savefig(to_path+'scatter'+str(n).zfill(4)+'.png',format='png',dpi=60)
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
               # exit()
        print "from path:", option.from_path
        print "to path:", option.to_path
        #if not os.path.exists(option.to_path):
        #        os.mkdir(option.to_path)
        main(option.from_path,option.to_path)

 
