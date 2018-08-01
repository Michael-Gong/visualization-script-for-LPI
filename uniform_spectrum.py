import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
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
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  


data = sdf.read('./Data/'+str(12).zfill(4)+".sdf",dict=True)
header=data['Header']
time=header['time']

px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
#field_ex = data['Particles/field_ex/subset_high_e/electron'].data/exunit
#field_ey = data['Particles/field_ey/subset_high_e/electron'].data/exunit
#field_bz = data['Particles/field_bz/subset_high_e/electron'].data/bxunit
gg = (px**2+py**2+1)**0.5

px = px [(abs(grid_y) < 3.2) & (gg > 0.0)]
py = py [(abs(grid_y) < 3.2) & (gg > 0.0)]
work_x = work_x [(abs(grid_y) < 3.2) & (gg > 0.0)]
work_y = work_y [(abs(grid_y) < 3.2) & (gg > 0.0)]

gg = (px**2+py**2+1)**0.5
theta = np.arctan2(py,px)*180.0/np.pi


x_1 = np.linspace(-12.,12.,200)
y_1,thrush = np.histogram(theta[work_x > work_y],bins=200,range=(-12.5,12.5))
y_2,thrush = np.histogram(theta[work_x < work_y],bins=200,range=(-12.5,12.5))


plt.subplot(2,1,1)
#    plt.subplot()
plt.scatter(theta[work_x > work_y], gg[work_x > work_y], c='magenta',  s=3, edgecolors='None', alpha=0.4, label='Work$_x$ > Work$_y$')
plt.scatter(theta[work_x < work_y], gg[work_x < work_y], c='cyan', s=3, edgecolors='None', alpha=0.6, label='Work$_x$ < Work$_y$')
#plt.legend(loc='upper left',fontsize=20,framealpha=0.5)

plt.xlim(-12,12)
#    plt.ylim(0,400)
plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
plt.ylabel('$\gamma$',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.ylim(25,800.0)
plt.twinx()

plt.plot(x_1,(y_1+y_2)/1000,'-k',linewidth=3,label='Total')
plt.plot(x_1,y_1/1000,linestyle='--',color='red',linewidth=3, label='Work$_x$ > Work$_y$')
plt.plot(x_1,y_2/1000,linestyle='--',color='blue',linewidth=3, label='Work$_x$ > Work$_y$')
plt.ylabel('dN/d'+r'$\theta$'+'[A.U.]', color='blue',fontdict=font)
plt.yticks(fontsize=20,color='blue')
plt.ylim(0,0.9)
plt.xlim(-12,12)


plt.subplot(2,1,2)

def pxpy_to_energy(gamma, weight):
    binsize = 200
    en_grid = np.linspace(0.5,799.5,200)
    en_bin  = np.linspace(0,800.0,201)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
      en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)

to_path='./figure_wrap_up/'

######### Parameter you should set ###########
start   =  12  # start time
stop    =  12  # end time
step    =  1  # the interval or step

#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
youwant =  ['electron_en']#,'electron_no_en']#,'ey','ex','ey_averaged','bz','bz_averaged','Subset_high_e_density','Subset_high_e_ekbar']
#youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
#youwant Derived electron_density,electron_ekbar...
#youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
if (os.path.isdir('jpg') == False):
  os.mkdir('jpg')
######### Script code drawing figure ################
for n in range(start,stop+step,step):
  #### header data ####
  data = sdf.read("./Data/"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  print('ok')
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  X, Y = np.meshgrid(x, y)
  
  for name in youwant:
      
      px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
      py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
      grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
      grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
      work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
      work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
      

      gg = (px**2+py**2+1.0)**0.5
      ww = data['Particles/Weight/subset_high_e/electron'].data*6.4e-6
      theta = np.arctan2(py,px)*180.0/np.pi
  

      px_x_d = px[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
      py_x_d = py[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
      gg_x_d = gg[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
      ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
      dist_x1, den1 = pxpy_to_energy(gg_x_d,ww_x_d)

      px_y_d = px[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
      py_y_d = py[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
      gg_y_d = gg[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
      ww_y_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x<=work_y)]
      dist_x2, den2 = pxpy_to_energy(gg_y_d,ww_y_d)

      px_d = px[ (grid_x>5) & (abs(grid_y)<3.) ]
      py_d = py[ (grid_x>5) & (abs(grid_y)<3.) ]
      gg_d = gg[ (grid_x>5) & (abs(grid_y)<3.) ]
      ww_d = ww[ (grid_x>5) & (abs(grid_y)<3.) ]
      dist_x, den = pxpy_to_energy(gg_d,ww_d)
      #plt.plot(dist_x*0.51,den,'-k',marker='o',markersize=10,linewidth=3,markevery=8)
      #plt.plot(dist_x1*0.51,den1,'-r',marker='^',markersize=10,linewidth=3,markevery=8)
      #plt.plot(dist_x2*0.51,den2,'-b',marker='s',markersize=10,linewidth=3,markevery=8)
      plt.plot(dist_x*0.51,den/2.0,'-k',linewidth=4,label='Total')
      plt.plot(dist_x1*0.51,den1/2.0,':r',linewidth=4,label='W$_x$ > W$_y$')
      plt.plot(dist_x2*0.51,den2/2.0,':b',linewidth=4,label='W$_x$ < W$_y$')

      print('g>200, LDA: ', np.sum(den1[dist_x1>200])/np.sum(den[dist_x>200]), ' TDA: ', np.sum(den2[dist_x2>200])/np.sum(den[dist_x>200]), 'charge:', np.sum(den[dist_x>200]))
      print('g>400, LDA: ', np.sum(den1[dist_x1>400])/np.sum(den[dist_x>400]), ' TDA: ', np.sum(den2[dist_x2>400])/np.sum(den[dist_x>400]), 'charge:', np.sum(den[dist_x>400]))
      print('g>600, LDA: ', np.sum(den1[dist_x1>600])/np.sum(den[dist_x>600]), ' TDA: ', np.sum(den2[dist_x2>600])/np.sum(den[dist_x>600]), 'charge:', np.sum(den[dist_x>600]))

      #### manifesting colorbar, changing label and axis properties ####
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.xlabel('Energy [MeV]',fontdict=font)
      plt.ylabel('dN/dE [MeV$^{-1}$]',fontdict=font)
      plt.xticks(fontsize=20); plt.yticks(fontsize=20);
      plt.yscale('log')
      plt.xlim(0,400)
      plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
      plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
      plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,
                wspace=None, hspace=None)


fig = plt.gcf()
fig.set_size_inches(9.2, 14)
fig.savefig('./figure_wrap_up/'+'comb_2_2_njp.png',format='png',dpi=160)
plt.close("all")

