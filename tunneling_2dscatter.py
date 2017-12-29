%matplotlib inline
#import sdf
import matplotlib
#matplotlib.use('agg')
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
font = {'family' : 'helvetica',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 12,
        }
part_number=20000
nsteps=2

insert1='./Dataqe/'
insert_n='_0'
#t=np.loadtxt(insert+'t'+'.txt')
z1=np.loadtxt(insert1+'z'+insert_n+'.txt')
y1=np.loadtxt(insert1+'y'+insert_n+'.txt')
x1=np.loadtxt(insert1+'x'+insert_n+'.txt')
#px=np.loadtxt(insert+'px'+'.txt')
#py=np.loadtxt(insert+'py'+'.txt')
#pz=np.loadtxt(insert+'pz'+'.txt')
#ey=np.loadtxt(insert+'e_part'+'.txt')
#bz=np.loadtxt(insert+'b_part'+'.txt')
#ay=np.loadtxt(insert+'a_part'+'.txt')
#radn=np.loadtxt(insert+'radn'+'.txt')
#radt=np.loadtxt(insert+'radt'+'.txt')
#opt=np.loadtxt(insert+'opt'+'.txt')
#eta=np.loadtxt(insert+'eta'+'.txt')

#t=np.reshape(t,(part_number,nsteps))
x1=np.reshape(x1,(part_number,nsteps))
y1=np.reshape(y1,(part_number,nsteps))
z1=np.reshape(z1,(part_number,nsteps))
#px=np.reshape(px,(part_number,nsteps))
#py=np.reshape(py,(part_number,nsteps))
#pz=np.reshape(pz,(part_number,nsteps))
#ey=np.reshape(ey,(part_number,nsteps))
#ay=np.reshape(ay,(part_number,nsteps))
#radn=np.reshape(radn,(part_number,nsteps))
#radt=np.reshape(radt,(part_number,nsteps))
#opt=np.reshape(opt,(part_number,nsteps))
#eta=np.reshape(eta,(part_number,nsteps))

#print(np.where(py[:,-1] > 0))

#gamma=np.sqrt(px**2+py**2+1)
insert1='./Datarr/'
insert_n='_0'
#t=np.loadtxt(insert+'t'+'.txt')
z2=np.loadtxt(insert1+'z'+insert_n+'.txt')
y2=np.loadtxt(insert1+'y'+insert_n+'.txt')
x2=np.loadtxt(insert1+'x'+insert_n+'.txt')

x2=np.reshape(x2,(part_number,nsteps))
y2=np.reshape(y2,(part_number,nsteps))
z2=np.reshape(z2,(part_number,nsteps))




makersize=0.8
plt.scatter(y1[:,-1]/2/np.pi,x1[:,-1]/2/np.pi,s=makersize,c=(192.0/255.0,0.0,0.0),label='QED',edgecolors='None')
plt.scatter(y2[:,-1]/2/np.pi,x2[:,-1]/2/np.pi,s=makersize,c=(0,192.0/255.0,0.0),label='LL',edgecolors='None')
#plt.scatter(x[:,-1]/2/np.pi,y[:,-1]/2/np.pi,s=makersize,c=(192.0/255.0,0.0,0.0),label='final',edgecolors='None')
#plt.legend(fontsize=14,loc='upper left',markerscale=5,scatterpoints=3,framealpha=0.0)
plt.xlim(-19.85,19.85)
plt.ylim(-0.95,5.95)
#ax.set_xlim([0.01,5.99])
#ax.set_ylim([-10.99,10.99])
#ax.set_zlim([-10.99,10.99])
ax = plt.gca() 
ax.invert_xaxis()
#ax.invert_yaxis()
plt.xlabel(r'Y\ [$\lambda_0$]',fontdict=font)
plt.ylabel(r'X\ [$\lambda_0$]',fontdict=font)
plt.xticks(fontsize=12); plt.yticks(fontsize=12);
#plt.title('electron scatter distribution',fontdict=font)
plt.grid(color='gray', linestyle='dotted', linewidth=0.5)
fig = plt.gcf()
fig.set_size_inches(5.4, 4.5)
fig.savefig('2dscatter_xy.png',format='png',dpi=640)
plt.close("all")
