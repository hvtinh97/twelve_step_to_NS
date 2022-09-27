# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:18:34 2022

@author: afdr9
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nx = 31
ny = 31
dx = 2 / (nx-1)
dy = 2 / (ny-1)
nt = 17
sigma = .25
nu = .05
dt = (sigma*dx*dy) / nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X, Y = np.meshgrid(x,y)

### Initial
u = np.ones((ny,nx))
u[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2

fig, axes = plt.subplots(3,2, figsize=(15,15), subplot_kw = {'projection':'3d'})
fig.suptitle('2D DIFFUSION', size = 20)

ax = axes[0][0]
ax.set_title('Initial', size = 15)
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.set_xlim(0,2)
ax.set_ylim(0,2)
ax.set_zlim(1,2)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)

### Diffusion equation
def diffuse(nt):
    for n in range(nt+1):
        un = np.copy(u)
        u[1:-1,1:-1] = un[1:-1,1:-1]\
                     + nu*(dt / dx**2)*(un[1:-1,2:] - 2*un[1:-1,1:-1] + un[1:-1,:-2])\
                     + nu*(dt / dy**2)*(un[2:,1:-1] - 2*un[1:-1,1:-1] + un[:-2,1:-1])
        u[0,:] = 1
        u[:,0] = 1
        u[-1,:] = 1
        u[:,-1] = 1
    ax.set_xlabel('x direction', size = 10)
    ax.set_ylabel('y direction', size = 10)
    ax.set_zlabel('z direction', size = 10)
    ax.set_xlim(0,2)
    ax.set_ylim(0,2)
    ax.set_zlim(1,2)
    ax.tick_params(axis = 'both', labelsize = 10)
    ax.plot_surface(X, Y, u, cmap = cm.viridis)                     

ax = axes[0][1]
diffuse(10)
ax.set_title('10 (s)', size = 15)

ax = axes[1][0]
diffuse(14)
ax.set_title('14 (s)', size = 15)

ax = axes[1][1]
diffuse(50)
ax.set_title('50 (s)', size = 15)

ax = axes[2][0]
diffuse(100)
ax.set_title('100 (s)', size = 15)

ax = axes[2][1]
diffuse(200)
ax.set_title('200 (s)', size = 15)





plt.show()


