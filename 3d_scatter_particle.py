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
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
        }
part_number=1000
nsteps=2

insert1='./Data/'
#t=np.loadtxt(insert+'t'+'.txt')
z1=np.loadtxt(insert1+'z'+'.txt')
y1=np.loadtxt(insert1+'y'+'.txt')
x1=np.loadtxt(insert1+'x'+'.txt')
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

insert1='./Data/'
#t=np.loadtxt(insert+'t'+'.txt')
z1=np.loadtxt(insert1+'z'+'.txt')
y1=np.loadtxt(insert1+'y'+'.txt')
x1=np.loadtxt(insert1+'x'+'.txt')
x1=np.reshape(x1,(part_number,nsteps))
y1=np.reshape(y1,(part_number,nsteps))
z1=np.reshape(z1,(part_number,nsteps))

insert2='./Datar_5/'
#t=np.loadtxt(insert+'t'+'.txt')
z2=np.loadtxt(insert2+'z'+'.txt')
y2=np.loadtxt(insert2+'y'+'.txt')
x2=np.loadtxt(insert2+'x'+'.txt')
x2=np.reshape(x2,(part_number,nsteps))
y2=np.reshape(y2,(part_number,nsteps))
z2=np.reshape(z2,(part_number,nsteps))


makersize=20
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x1[:,-1]/2/np.pi,y1[:,-1]/2/np.pi,z1[:,-1]/2/np.pi,s=makersize,depthshade=True,c=(192/255.0,0,0),edgecolor='none')
#ax.scatter(x1[:,0]/2/np.pi,y1[:,0]/2/np.pi,z1[:,0]/2/np.pi,s=makersize,depthshade=True,c=(192/255.0,0,0),edgecolor='none')
ax.scatter(x2[:,-1]/2/np.pi,y2[:,-1]/2/np.pi,z2[:,-1]/2/np.pi,s=makersize,depthshade=True,c=(0,192/255.0,0),edgecolor='none')
#ax.scatter(x2[:,0]/2/np.pi,y2[:,0]/2/np.pi,z2[:,0]/2/np.pi,s=makersize,depthshade=True,c=(192/255.0,0,0),edgecolor='none')



ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_xlim([0,100])
ax.set_ylim([-50,50])
ax.set_zlim([-50,50])
ax.view_init(elev=30, azim=-60)


