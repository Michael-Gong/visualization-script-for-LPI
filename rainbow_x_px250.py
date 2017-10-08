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
        'size'   : 20,  
        }  



######### Parameter you should set ###########
start   =  1  # start time
stop    =  400  # end time
step    =  1  # the interval or step

Data = 'Data0/'
name = 'electron scattering plot'
amplitude = '250'
######### Script code drawing figure ################

data0 = sdf.read("./Dataa"+amplitude+"no/0000.sdf",dict=True)
#data1 = sdf.read("./Dataa350rr/"+str(n).zfill(4)+".sdf",dict=True)
header=data0['Header']
time=header['time']
x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
id0 = data0['Particles/ID/electron'].data
pos0= np.zeros(5000)
print(np.max(id0-1),np.min(id0-1),x0[id0-1])
for i in range(5000):
    pos0[id0[i]-1]=y0[i]

data0 = sdf.read("./Dataa"+amplitude+"rr/0000.sdf",dict=True)
#data1 = sdf.read("./Dataa350rr/"+str(n).zfill(4)+".sdf",dict=True)
header=data0['Header']
time=header['time']
x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
id0 = data0['Particles/ID/electron'].data
pos1= np.zeros(5000)
print(np.max(id0-1),np.min(id0-1),x0[id0-1])
for i in range(5000):
    pos1[id0[i]-1]=y0[i]

data0 = sdf.read("./Dataa"+amplitude+"qe/0000.sdf",dict=True)
#data1 = sdf.read("./Dataa350rr/"+str(n).zfill(4)+".sdf",dict=True)
header=data0['Header']
time=header['time']
x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
id0 = data0['Particles/ID/electron'].data
pos2= np.zeros(5000)
print(np.max(id0-1),np.min(id0-1),x0[id0-1])
for i in range(5000):
    pos2[id0[i]-1]=y0[i]
    
def main(from_path, to_path):
  for n in range(start,stop+step,step):
    if n < 40:
        x_min,x_max=45,55
        y_min,y_max=500,1100
    elif n < 50:
        x_min,x_max=45,55
        y_min,y_max=0,1100
    elif n < 70: 
        x_min,x_max=45,55
        y_min,y_max=-600,1100
    elif n < 90: 
        x_min,x_max=40,60
        y_min,y_max=-1000,1100
    elif n < 125: 
        x_min,x_max=40,60
        y_min,y_max=-1000,1100
    elif n < 180: 
        x_min,x_max=35,65
        y_min,y_max=-1000,1100
    elif n < 235: 
        x_min,x_max=30,70
        y_min,y_max=-400,400
    elif n < 290: 
        x_min,x_max=25,75
        y_min,y_max=-400,400
    else:
        x_min,x_max=25-(n*1.0-290.0)*0.1,75+(n*1.0-290.0)*0.1
        y_min,y_max=-300,300
    #plt.xlim(x_min,x_max)
    #plt.ylim(y_min,y_max)

    datalaser = sdf.read("./Datalaser/"+str(n).zfill(4)+".sdf",dict=True)
    laser_x   = datalaser['Grid/Grid_mid'].data[0]/1.0e-6
    laser_y   = datalaser['Grid/Grid_mid'].data[1]/1.0e-6
    X,Y = np.meshgrid(laser_x,laser_y)
    ex  = datalaser['Electric Field/Ey'].data/exunit/250.0
    if np.min(ex.T) == np.max(ex.T):
          continue
    eee=np.max([-np.min(ex.T),np.max(ex.T)])
    levels = np.linspace(-eee, eee, 24)
    #### header data ####
    plt.subplots_adjust(left=0.05,right=0.9,bottom=0.1,top=0.95,wspace=0.15,hspace=0.2)

    plt.subplot(1,3,1)
    data0 = sdf.read("./Dataa"+amplitude+"no/"+str(n).zfill(4)+".sdf",dict=True)
    #data1 = sdf.read("./Dataa350rr/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    px0 = data0['Particles/Px/electron'].data/(m0*v0)
    id0 = data0['Particles/ID/electron'].data

#    x   = data0['Grid/Grid_mid'].data[0]/1.0e-6
#    y   = data0['Grid/Grid_mid'].data[1]/1.0e-6
#  X,Y = np.meshgrid(x,y)
#  ex  = data0['Electric Field/Ey'].data/exunit
#  if np.min(ex.T) == np.max(ex.T):
#        continue
#  eee=np.max([-np.min(ex.T),np.max(ex.T)])
#  levels = np.linspace(-eee, eee, 24)
#  plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.RdGy)
    plt.scatter(x0,px0,s=20,c=abs(pos0[id0-1]),cmap=cm.rainbow,label='No RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=2,fontsize=20.0)
    plt.xlim(40,90)
    plt.ylim(-50,1100)
    #plt.text(x_min+0.1*(x_max-x_min),y_max-0.1*(y_max-y_min),r'$\xi_0$='+amplitude+'   No RR',fontsize=25)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Px [$m_ec$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=12);
    plt.title(name+' at '+str(round(time/3.333e-15,2))+' $T_0$',fontdict=font)

    plt.subplot(1,3,2)
    data0 = sdf.read("./Dataa"+amplitude+"rr/"+str(n).zfill(4)+".sdf",dict=True)
#    data1 = sdf.read("./epoch2drr/Data1/"+str(n).zfill(4)+".sdf",dict=True)
#    data2 = sdf.read("./epoch2dqe/Data1/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    px0 = data0['Particles/Px/electron'].data/(m0*v0)
    id0 = data0['Particles/ID/electron'].data
#    plt.scatter(x1,y1,s=8,c=(0,192.0/255.0,0),label='LL RR',edgecolors='None') 
    #plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.bwr)
    plt.scatter(x0,px0,s=20,c=abs(pos1[id0-1]),cmap=cm.rainbow,label='LL RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=2,fontsize=20.0)
    plt.xlim(40,90)
    plt.ylim(-50,1100)
    #plt.text(x_min+0.1*(x_max-x_min),y_max-0.1*(y_max-y_min),r'$\xi_0$='+amplitude+'   LL RR',fontsize=25)    
    #plt.text(5,45,r'$\xi_0=350$',fontsize=20)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Px [$m_ec$]',fontdict=font)    
    plt.xticks(fontsize=20); plt.yticks(fontsize=12);
    plt.title(name+' at '+str(round(time/3.333e-15,2))+' $T_0$',fontdict=font)

    plt.subplot(1,3,3)
    data0 = sdf.read("./Dataa"+amplitude+"qe/"+str(n).zfill(4)+".sdf",dict=True)
    header=data0['Header']
    time=header['time']
    x0  = data0['Grid/Particles/electron'].data[0]/1.0e-6
    y0  = data0['Grid/Particles/electron'].data[1]/1.0e-6
    px0 = data0['Particles/Px/electron'].data/(m0*v0)
    id0 = data0['Particles/ID/electron'].data
#    plt.scatter(x1,y1,s=8,c=(0,192.0/255.0,0),label='LL RR',edgecolors='None') 
    #plt.contourf(X, Y, ex.T, levels=levels, cmap=cm.bwr)
    plt.scatter(x0,px0,s=20,c=abs(pos2[id0-1]),cmap=cm.rainbow,label='QED RR',edgecolors='None') 
    plt.legend(loc='upper right',framealpha=1.0,markerscale=2,fontsize=20.0)
    #plt.text(x_min+0.1*(x_max-x_min),y_max-0.1*(y_max-y_min),r'$\xi_0$='+amplitude+'   QED RR',fontsize=25)
    plt.xlim(40,90)
    plt.ylim(-50,1100)
    #plt.text(5,45,r'$\xi_0=350$',fontsize=20)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Px [$m_ec$]',fontdict=font)    
    plt.xticks(fontsize=20); plt.yticks(fontsize=12);
    plt.title(name+' at '+str(round(time/3.333e-15,2))+' $T_0$',fontdict=font)
    
    
    #plt.subplots_adjust(left=0.05,right=0.85,bottom=0.1,top=0.95,wspace=0.15,hspace=0.2)
    cbar=plt.colorbar(ticks=[0.5, 1.0, 1.5, 2.0],cax=plt.axes([0.94,0.1,0.01,0.85]))
    cbar.set_label('Initial transverse position $y_0$',fontdict=font) 
    fig = plt.gcf()
    fig.set_size_inches(33, 9)
    fig.savefig('./jpg_x_px250/gif'+str(n).zfill(4)+'.png',format='png',dpi=45)
    plt.close("all")
    fig.show()
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
