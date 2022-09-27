# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 18:03:38 2022

@author: afdr9
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nx = 50
ny = 50
dx = 2 / (nx-1)
dy = 1 / (ny-1)
nt = 100

x = np.linspace(0,2,nx)
y = np.linspace(0,1,ny)
X, Y = np.meshgrid(x,y)
b = np.zeros((ny,nx))
b[int(ny/4), int(nx/4)] = 100
b[int(3*ny/4), int(3*nx/4)] = -100
p = np.zeros((ny,nx))
pd = np.zeros((ny,nx))

for i in range(nt):
    pd = np.copy(p)
    p[1:-1,1:-1] = ((pd[1:-1,2:] + pd[1:-1,0:-2])*(dy**2)\
                 + (pd[2:,1:-1] + pd[0:-2,1:-1])*(dx**2)\
                 - (b[1:-1,1:-1]*(dx**2)*(dy**2)))\
                 / (2*(dx**2 + dy**2))
    p[0,:] = 0
    p[-1,:] = 0
    p[:,0] = 0
    p[:,-1] = 0

def plot3d(x, y, p):
    X, Y = np.meshgrid(x,y)
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x direction', size = 10)
    ax.set_ylabel('y direction', size = 10)
    ax.set_zlabel('z direction', size = 10)
    ax.tick_params(axis = 'both', size = 10)
    ax.view_init(30,225)
    ax.plot_surface(X, Y, p, cmap = cm.viridis)
plot3d(x,y,p)
plt.show()
