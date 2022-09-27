# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 22:47:21 2022

@author: afdr9
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nx = 41
ny = 41
dx = 2 / (nx-1)
dy = 2 / (ny-1)
nt = 120
sigma = .0009
nu = 0.01
dt = (sigma*dx*dy) / nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X, Y = np.meshgrid(x,y)

### Initial
u = np.ones((ny,nx))
u[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2
v = np.ones((ny,nx))
v[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2

fig, axes = plt.subplots(2,2, figsize=(15,15), subplot_kw = {'projection':'3d'})
fig.suptitle("2D BURGERS' EQUATION", size = 20)

ax = axes[0][0]
ax.set_title('Initial values for u', size = 15)
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)

ax = axes[0][1]
ax.set_title('Initial values for v', size = 15)
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, v, cmap = cm.viridis)

### u field & v field
for n in range(nt+1):
    un = np.copy(u)
    vn = np.copy(u)
    
    u[1:-1,1:-1] = un[1:-1,1:-1]\
             - (dt/dx)*un[1:-1,1:-1]*(un[1:-1,1:-1] - un[1:-1,:-2])\
             - (dt/dy)*vn[1:-1,1:-1]*(un[1:-1,1:-1] - un[:-2,1:-1])\
             - (nu*dt / dx**2)*(un[1:-1,2:] - 2*un[1:-1,1:-1] + un[1:-1,:-2])\
             - (nu*dt / dy**2)*(un[2:,1:-1] - 2*un[1:-1,1:-1] + un[:-2,1:-1])
    v[1:-1,1:-1] = vn[1:-1,1:-1]\
             - (dt/dx)*vn[1:-1,1:-1]*(vn[1:-1,1:-1] - vn[1:-1,:-2])\
             - (dt/dy)*vn[1:-1,1:-1]*(vn[1:-1,1:-1] - vn[:-2,1:-1])\
             - (nu*dt / dx**2)*(vn[1:-1,2:] - 2*vn[1:-1,1:-1] + vn[1:-1,:-2])\
             - (nu*dt / dy**2)*(vn[2:,1:-1] - 2*vn[1:-1,1:-1] + vn[:-2,1:-1])

    u[0,:] = 1
    u[:,0] = 1
    u[-1,:] = 1
    u[:,-1] = 1
    v[0,:] = 1
    v[:,0] = 1
    v[-1,:] = 1
    v[:,-1] = 1
    
ax = axes[1][0]
ax.set_title('u field', size = 15)
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)

ax = axes[1][1]
ax.set_title('u field', size = 15)
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, v, cmap = cm.viridis)
plt.show()
