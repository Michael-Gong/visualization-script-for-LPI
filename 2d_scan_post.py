#%matplotlib inline
import sdf
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
import sys
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
       }

font_size = 25

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.viridis(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_viridis = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.rainbow(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_rainbow = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def pxpy_to_energy(gamma, weight):
    binsize = 200
    en_grid = np.linspace(50,19950,200)
    en_bin  = np.linspace(0,20000.0,201)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)


def theta_to_grid(theta, weight):
    binsize = 240
    theta_grid = np.linspace(-119.5,119.5,240)
    theta_bin  = np.linspace(-120,120,241)
    theta_value = np.zeros_like(theta_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            theta_value[i] = sum(weight[ (theta_bin[i]<=theta) & (theta<theta_bin[i+1]) ])
    return (theta_grid, theta_value)


def cal_ionization_status(from_dir,nth):
    data = sdf.read('./'+from_dir+'/tot'+str(nth).zfill(4)+'.sdf',dict=True)
    header=data['Header']
    time=header['time']
    print('processing '+from_dir+'; time:'+str(nth))
    #for i in range(np.size(youwant)):
    if 'Particles/ID/subset_sub/Ion1' in data:
        if 'Particles/ID/subset_sub/Ion' in data:
          num_C  = np.sum(data['Particles/ID/subset_sub/Ion'].data*data['Particles/Weight/subset_sub/Ion'].data)
          num_C1 = np.sum(data['Particles/ID/subset_sub/Ion1'].data*data['Particles/Weight/subset_sub/Ion1'].data)
          return (num_C1)/(num_C1+num_C)
        else:
          return 1.0
    else:
        return 0.0
    
def plot_2d_scan(x,y,data,time,n):
    X, Y = np.meshgrid(x, y, indexing='ij')
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    levels = np.linspace(0,100, 101)
    img = ax.pcolormesh(X, Y, data, norm=colors.Normalize(vmin=0, vmax=100), cmap=mycolor_jet,edgecolors='k')
#  cax = fig.add_axes([0.65,0.94,0.25,0.02])
#  cbar=fig.colorbar(img,cax=cax, ticks=[1e2,1e5],orientation='horizontal')
    cbar=fig.colorbar(img,pad=0.01, ticks=[0,50,100])
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    cbar.set_label('ionized fraction [%]',fontdict=font)
      #ax.set_xlim(10,50)
      #ax.set_ylim(0.,1.)
    ax.set_xlabel(r'log$_{10}$$a_0$',fontdict=font)
    ax.set_ylabel(r'log$_{10}$$\varepsilon$[eV]',fontdict=font)
#  ax.set_yscale('log')
    ax.set_ylim(0,5)
    ax.set_xlim(-2,3)
  #ax.set_xticklabels([0,90,180,270])
  #ax.set_yticklabels([0.1,1,10,100,1000])

  #ax.set_theta_zero_location('N')
    #  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
    ax.tick_params(axis='x',labelsize=font_size) 
    ax.tick_params(axis='y',labelsize=font_size)
    ax.set_title('time t = '+str(time)+' fs', fontsize=20)
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)

#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
#  plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
#  plt.ylabel('time [fs]',fontdict=font)
#  plt.xticks([-135,-90,-45,0,45,90,135],fontsize=font_size); plt.yticks([0,500,1000,1500],fontsize=font_size);
#  plt.title(r'$dN/d\theta$'+' for no RR', fontsize=font_size)
#  plt.xlim(-120,120)
#  plt.ylim(0,1650)
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    plt.subplots_adjust(top=0.96, bottom=0.13, left=0.16, right=0.95, hspace=0.10, wspace=0.05)

    fig = plt.gcf()
    fig.set_size_inches(10., 8.5)
#fig.set_size_inches(5, 4.5)
    fig.savefig('./2d_scan_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    return 0


if __name__ == "__main__":
  grid_a0 = np.linspace(1,21,21)
  grid_e0 = np.linspace(1,21,21)
  from_path = './'
  start   = 0
  stop    = 100 
  step    = 1
  grid_t0 = np.linspace(start,stop,(stop-start+step)/step)*10.0

  grid_data = np.zeros([21,21,(stop-start+step)/step])
  for n in range(start,stop+step,step):
      for i in range(np.size(grid_a0)):
          for j in range(np.size(grid_e0)):
              path_dir = 'Dataf'+str(int(grid_a0[i]))+'e'+str(int(grid_e0[j]))
              grid_data[i,j,(n-start)/step] = cal_ionization_status(path_dir,n)*100.0
              print(path_dir+' fraction:',grid_data[i,j,(n-start)/step])
   
      plot_2d_scan(-2+(grid_a0-1)/4,(grid_e0-1.0)/4,grid_data[:,:,(n-start)/step],grid_t0[(n-start)/step],n)
      print('finish 2d_scan figure for n='+str(n))
          
          
