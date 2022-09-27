# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 22:57:33 2022

@author: afdr9
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nx = 41
ny = 41
nt = 500
dt = .001
nu = .1

p = np.zeros((ny,nx))
b = np.zeros((ny,nx))
rho = 1
u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
dx = 2 / (nx-1)
dy = 2 / (ny-1)

def build_up_b(b, rho, dt, u, v, dx, dy):
    b[1:-1,1:-1] = ((rho/dt)*(((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)) + ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))))\
                 - (((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)))\
                 - (2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy))*((v[1:-1,2:]-v[1:-1,0:-2])/(2*dx)))\
                 - (((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))*((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy)))
    return b

nit = 50
def poisson(p,dx,dy,b):
    for i in range(nit):
        pn = np.copy(p)       
        p[1:-1,1:-1] = ((((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2)\
                     + (pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)\
                     / (2*(dx**2 + dy**2)))\
                     - ((rho * dx**2 * dy**2)/(2*(dx**2 + dy**2)))\
                     * b[1:-1,1:-1]
        p[:,-1] = p[:,-2]
        p[:,0]  = p[:,1]
        p[0,:]  = p[1,:]
        p[-1,:] = 0
    return p

def cavity(u,v,nt,dt,dx,dy,p,rho,nu):
    b = np.zeros((ny,nx))
    for n in range(nt):
        un = np.copy(u)
        vn = np.copy(v)
        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = poisson(p,dx,dy,b)
        u[1:-1,1:-1] = (un[1:-1,1:-1])\
                     - (un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1] - un[1:-1,0:-2]))\
                     - (vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1]-un[0:-2, 1:-1]))\
                     - ((dt/(2*rho*dx))*(p[1:-1,2:]-p[1:-1,0:-2]))\
                     + (nu*((dt/(dx**2))*(un[1:-1,2:] - 2*un[1:-1,1:-1] + un[1:-1,0:-2])\
                     + (dt/(dy**2))*(un[2:,1:-1] - 2*un[1:-1,1:-1] + un[0:-2,1:-1])))
        v[1:-1,1:-1] = (vn[1:-1,1:-1])\
                     - (un[1:-1,1:-1]*(dt/dx)*(vn[1:-1,1:-1] - vn[1:-1,0:-2]))\
                     - (vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1]-un[0:-2, 1:-1]))\
                     - ((dt/(2*rho*dy))*(p[2:,1:-1]-p[0:-2,1:-1]))\
                     + (nu*((dt/(dx**2))*(vn[1:-1,2:] - 2*vn[1:-1,1:-1] + vn[1:-1,0:-2])\
                     + (dt/(dy**2))*(vn[2:,1:-1] - 2*vn[1:-1,1:-1] + vn[0:-2,1:-1])))                     
        u[0,:]  = 0
        u[:,0]  = 0
        u[:,-1] = 0
        u[-1,:] = 2
        v[0,:]  = 0
        v[:,0]  = 0
        v[:,-1] = 0
        v[-1,:] = 0
    return u, v, p        
        
nt = 100
u, v, p = cavity(u,v,nt,dt,dx,dy,p,rho,nu)

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X, Y = np.meshgrid(x,y)

plt.figure(figsize=(10,10))
plt.contourf(X, Y, p, level = 30, cmap ='jet')
plt.colorbar()
plt.contour(X, Y, p, level = 30, cmap = 'jet')
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
plt.streamplot(X, Y, u, v, cmap = 'jet')
plt.show()




      
        
