# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:50:56 2022

@author: afdr9
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nx = 81
ny = 81
dx = 2 / (nx-1)
dy = 2 / (ny-1)
nt = 100
sigma = .2
dt = sigma*dx
c = 1

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X, Y = np.meshgrid(x,y)

### Initial
u = np.ones((ny,nx))
u[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2

fig, axes = plt.subplots(2,2, figsize = (15,15), subplot_kw = {'projection':'3d'})
fig.suptitle('2D LINEAR CONVECTION', size = 20)

ax = axes[0][0]
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.set_title('Initial', size = 15)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)

### Method 1
u = np.ones((ny,nx))
u[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2
for n in range(nt+1):
    un = np.copy(u)
    row, col = np.shape(u)
    for j in range(row):
        for i in range(col):
            u[j,i] = un[j,i] - c*(dt/dx)*(un[j,i] - un[j,i-1])\
                   - c*(dt/dy)*(un[j,i] - un[j-1,i])
            u[0,:] = 1
            u[:,0] = 1
            u[-1,:] = 1
            u[:,-1] = 1

ax = axes[0][1]
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.set_title('Method 1', size = 15)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)

### Method 2
u = np.ones((ny,nx))
u[int(.5/dy):int(1/dy + 1), int(.5/dx):int(1/dx + 1)] = 2
for n in range(nt+1):
    u[1:,1:] = un[1:,1:] - c*(dt/dx)*(un[1:,1:] - un[1:,:-1])\
                   - c*(dt/dy)*(un[1:,1:] - un[:-1,1])
    u[0,:] = 1
    u[:,0] = 1
    u[-1,:] = 1
    u[:,-1] = 1

ax = axes[1][0]
ax.set_xlabel('x direction', size = 10)
ax.set_ylabel('y direction', size = 10)
ax.set_zlabel('z direction', size = 10)
ax.set_title('Method 2', size = 15)
ax.tick_params(axis = 'both', labelsize = 10)
ax.plot_surface(X, Y, u, cmap = cm.viridis)
plt.show()





