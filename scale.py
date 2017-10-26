import cmath
import numpy as np  
from matplotlib import pyplot as plt  
from matplotlib import animation 

# This is the relation between dz, z_max, L 
def relation_dz(L, z_max, ns):
    return L*np.log(ns/(2-1.0/np.sqrt(1.0+0.064**2)*(2.0-ns*np.exp(-z_max/L))))

def cal_main_solve(ns, L, dz): #this is analytical method. Finally fail!!!
    a = 2.0
    b = ns
    c = L
    d = (2.0-ns*np.exp(0.0-dz/L))**2
    z = dz
#    print(a,b,c,d,z)
    const = -1j*c*np.sqrt(d)*np.log(2.0*np.exp(z/c)*(-np.sqrt(0j+(a**2-2*a*b*np.exp(-z/c)+b**2*np.exp(-2*z/c)-d)/d)+1j*(-a**2+a*b*np.exp(-z/c)+d)/np.sqrt(d*(d-a**2)+0j)))/np.sqrt(d-a**2+0j)

    z = z_max
    return 1j*c*np.sqrt(d)*np.log(2.0*np.exp(z/c)*(-np.sqrt(0j+(a**2-2*a*b*np.exp(-z/c)+b**2*np.exp(-2*z/c)-d)/d)+1j*(-a**2+a*b*np.exp(-z/c)+d)/np.sqrt(d*(d-a**2)+0j)))/np.sqrt(d-a**2+0j)+const

def deviation(y_f, l_y):
    return y_f-l_y

def cal_main_solve_2(ns, L, dz, ly):
    dd = (ly-0.0)/10000
    y = 0.0
    z = dz 
    for y in np.arange(10000):
        dy1 = dd
        dz1 = dy1*np.sqrt((2.0-ns*np.exp(-z/L))**2/(2.0-ns*np.exp(-dz/L))**2/np.sin(np.pi/2.0/90.0*89.99)**2-1.0)
        y = y + dy1
        z = z + dz1   
    return z  


if __name__ == '__main__':
    l_y = 7.0
    z_max = 7.3
    ns = 100.0
    for L in np.linspace(0.01,2,500):
        dz = relation_dz(L, z_max, ns)
 #       print(L, dz)
        z = cal_main_solve_2(ns, L, dz, l_y)
        #print("L=",L," and deviation delta=",deviation(y_f, l_y))
        print("L=",L," and deviation delta=",z-z_max)
        #print('L=',L,' and deviation delta=')

