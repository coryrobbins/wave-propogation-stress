# 3D wave propogation stress modeling
# Differing materials / fluids during propogation
# Python version of MatLab code - converted for a close colleague

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


''' Custom settings - number of cells nn, tt '''
nn = 500
tt = 600        


''' Initialize matricies '''
v = np.mat(np.zeros(shape = (nn, tt)))
e = np.mat(np.zeros(shape = (nn, tt)))
Vp = np.mat(np.zeros(shape = (nn, tt)))
Vm = np.mat(np.zeros(shape = (nn, tt)))
r = np.mat(np.zeros(shape = (nn, tt)))
cp = np.mat(np.zeros(shape = (nn, tt)))
stress = np.mat(np.zeros(shape = (nn, tt)))
Sm = np.mat(np.zeros(shape = (nn, tt)))
Sp = np.mat(np.zeros(shape = (nn, tt)))

tt = (tt-1)   # Because Python uses 0 as first row in matrices
nn = (nn-1)

e[99:149, 0] = 2
r[0:199, 0:tt]= 1       # Density matrix
r[199:nn, 0:tt] = 0.8
cp[0:199, 0:tt] = 1     # Speed
cp[199:nn, 0:tt] = 0.8 
alpha = 1


''' Begin Calculations '''
for t in range(0, tt-1):  
    
    ''' Boundary Conditions '''
    v[nn-1, 0:tt] = 0
    e[nn, 0:tt] = (e[nn-1, tt] + r[nn, tt] * cp[nn, tt] * ( 2 * v[nn-1, tt ] - v[nn-2, tt] ))
    
    
    for n in range(1, nn-2):
        
        '''V Minus '''
        Vm[n,t] = ((r[n,t]*cp[n,t]*cp[n,t]*e[n,t]-r[n-1,t]*cp[n-1,t]*cp[n-1,t]*e[n-1,t])-r[n-1,t]*
                   cp[n-1,t]*(v[n,t]-v[n-1,t])) /(r[n-1,t]*cp[n-1,t]+r[n,t]*cp[n,t])
        
        ''' V Plus '''      
        Vp[n,t]= ((r[n+1,t]*cp[n+1,t]*cp[n+1,t]*e[n+1,t] - r[n,t]*cp[n,t]*cp[n,t]*e[n,t])+r[n+1,t]*
                  cp[n+1,t]*(v[n+1,t]- v[n,t])) / (r[n,t]*cp[n,t]+r[n+1,t]*cp[n+1,t])

    for n in range(0, nn-2):

        Sp[n,t] = (r[n,t] * cp[n,t] * Vp[n,t])
        ''' Sigma Minus '''

        Sm[n,t]= -(r[n,t]*cp[n,t]*Vm[n,t]) 
        ''' Sigma Plus '''
        
    for n in range(0, nn-2):    
        
        ''' Pv equation to find velocity '''
        v[n,t+1] = (v[n,t]) + alpha *  (Sp[n,t]-Sm[n,t])/ r[n,t]
        e[n,t+1] = e[n,t] + alpha * (Vp[n,t] - Vm[n,t])
        
        ''' Final stress calculation '''
        stress[n,t] = (e[n,t] * r[n,t] * cp[n,t] * cp[n,t]) 
       
        
''' Begin Plotting Section '''
fig = plt.figure()
ax = fig.gca(projection='3d')

# Set up mesh with X,Y = 500x500 grid
X = np.arange(0, nn, 1)
Y = np.arange(0, tt, 1)
X, Y = np.meshgrid(X, Y)   # mesh together

Z = (stress[X, Y])

# Now plot surface
ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

# Display plot
plt.show()
