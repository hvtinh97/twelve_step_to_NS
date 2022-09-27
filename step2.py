# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 20:25:42 2022

@author: afdr9
"""

#import numpy as np
#import matplotlib.pyplot as plt
#
#nx = 41
#dx = 2 / (nx-1)
#nt = 25
#dt = .025
#
#u = np.ones(nx)
#u[int(.5/dx) : int(1/dx + 1)] = 2
#plt.style.use(['bmh'])
#plt.figure(figsize = (10,10))
#plt.title('Step 2', size = 15)
#plt.xlabel('Distance', size = 15)
#plt.ylabel('Velocity', size = 15)
#plt.tick_params(axis = 'both', labelsize = 15)
#plt.plot(np.linspace(0,2,nx), u, 'o:', color = 'blue', lw = 1.5, ms = 6, label = 'Initial value')
#
#un = np.ones(nx)
#for n in range(nt):
#    un = u.copy()
#    for i in range(1,nx):
#        u[i] = un[i] - un[i] * (dt/dx) * (un[i] - un[i-1])
#plt.plot(np.linspace(0,2,nx), u, 'o:', color = 'red', lw = 1.5, ms = 6, label = 'Solution')
#plt.legend(loc = 'upper right', fontsize = 15)
#plt.show()

###Apply function
#import numpy as np
#import matplotlib.pyplot as plt
#
#plt.style.use(['bmh'])
#fig, axes = plt.subplots(2,2, figsize = (10,10))
#fig.suptitle('APPLY FUNCTION', size = 20)
#
#def linconv(nx):
#    dx = 2 / (nx-1)
#    nt = 20
#    dt = .025
#    c = 1
#    
#    u = np.ones(nx)
#    u[int(.5/dx) : int(1/dx + 1)] = 2
#    
#    un = np.ones(nx)
#    for n in range(nt):
#        un = u.copy()
#        for i in range(nx):
#            u[i] = un[i] - c * (dt/dx) * (un[i] - un[i-1])
#    ax.set_xlabel('Distance', size = 10)
#    ax.set_ylabel('Velocity', size = 10)
#    ax.set_xlim(0,2)
#    ax.set_ylim(-1,4)
#    ax.tick_params(axis = 'both', labelsize = 10)
#    ax.plot(np.linspace(0,2,nx), u, 'o:', color = 'blue', lw = 1.5, ms = 6)
#
#ax = axes[0][0]
#linconv(41)
#ax.set_title('41 nodes', size = 15)
#
#ax = axes[0][1]
#linconv(61)
#ax.set_title('41 nodes', size = 15)
#
#ax = axes[1][0]
#linconv(61)
#ax.set_title('41 nodes', size = 15)
#    
#ax = axes[1][1]
#linconv(85)
#ax.set_title('41 nodes', size = 15) 
#plt.show()
    
###CFL condition
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['bmh'])
fig, axes = plt.subplots(3,2, figsize = (10,20))
fig.suptitle('CFL CONDITION', size = 20)

def linconv(nx):
    dx = 2 / (nx-1)
    nt = 20
    c = 1
    sigma = .5
    dt = sigma*dx
    
    u = np.ones(nx)
    u[int(.5/dx) : int(1/dx + 1)] = 2
    
    un = np.ones(nx)
    for n in range(nt):
        un = u.copy()
        for i in range(nx):
            u[i] = un[i] - c * (dt/dx) * (un[i] - un[i-1])
    ax.set_xlabel('Distance', size = 10)
    ax.set_ylabel('Velocity', size = 10)
    ax.set_xlim(0,2)
    ax.set_ylim(0.5,2.5)
    ax.tick_params(axis = 'both', labelsize = 10)
    ax.plot(np.linspace(0,2,nx), u, 'o:', color = 'red', lw = 1.5, ms = 6)

ax = axes[0][0]
linconv(41)
ax.set_title('41 nodes', size = 15)

ax = axes[0][1]
linconv(61)
ax.set_title('61 nodes', size = 15)

ax = axes[1][0]
linconv(81)
ax.set_title('81 nodes', size = 15)

ax = axes[1][1]
linconv(101)
ax.set_title('101 nodes', size = 15)

ax = axes[2][0]
linconv(121)
ax.set_title('121 nodes', size = 15)

plt.show()










    
    
    
    
    