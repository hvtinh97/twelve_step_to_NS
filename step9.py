# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 19:50:00 2022

@author: afdr9
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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
    
def laplace2d(l1norm_target):
    l1norm = 1
    while l1norm > l1norm_target:
        pn = np.copy(p)
        p[1:-1,1:-1] = (((dy**2)*(pn[1:-1, 2:] + pn[1:-1, 0:-2]))\
                     + ((dx**2)*(pn[2:,1:-1] + pn[0:-2,1:-1])))\
                     / (2*(dx**2 + dy**2))
        p[:,0] = 0
        p[:,-1] = y
        p[0,:] = p[1,:]
        p[-1,:] = p[-2,:]
        l1norm = (np.sum(np.abs(p[:]) - np.abs(pn[:]))) / (np.sum(np.abs(pn[:])))
    
nx = 31
ny = 31
dx = 2 / (nx-1)
dy = 1 / (nx-1)
c = 1

x = np.linspace(0,2,nx)
y = np.linspace(0,1,ny)

p = np.zeros((ny,nx))
p[:,0] = 0
p[:,-1] = y
p[0,:] = p[1,:]
p[-1,:] = p[-2,:]

plot3d(x, y, p)
laplace2d(1e-4)
plot3d(x, y, p)
plt.show()
