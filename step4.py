# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:25:22 2022

@author: afdr9
"""
#import numpy as np
#import sympy as smp
#import matplotlib.pyplot as plt
#from sympy import init_printing
#init_printing(use_latex = True, forecolor = 'Yellow')
#
#x, nu, t = smp.symbols('x  nu t', real = True)
#phi = smp.exp((-(x - 4*t)**2) / (4*nu*(t+1)))\
#    + smp.exp((-(x- 4*t -2*smp.pi)**2) / (4*nu*(t+1)))
#print('phi =')
#display(phi)
#
#phiprime = smp.diff(phi,x)
#print('*****\nphiprime =')
#display(phiprime)
#
#u = -((2*nu)/phi)*phiprime + 4
#print('*****\nu =')
#display(u)
#
#from sympy.utilities.lambdify import lambdify
#ufunc = lambdify((nu, t, x), u)
#print('*****\nufunc = ', ufunc(3,1,4))
#
####Initial value
#nx = 101
#dx = (2*np.pi) / (nx-1)
#nt = 100
#nu = .07
#dt = nu*dx
#t = 0
#
#x = np.linspace(0, 2*np.pi, nx)
#u = np.asarray([ufunc(nu, t, x0)for x0 in x])
#print('*****\nu =', u)
#
#plt.style.use(['bmh'])
#plt.figure(figsize = (10,10))
#plt.title("BURGER'S EQUATION", size = 20)
#plt.xlabel('Distance', size = 10)
#plt.ylabel('Velocity', size = 10)
#plt.xlim(0,7)
#plt.ylim(0,8)
#plt.tick_params(axis = 'both', labelsize = 10)
#plt.plot(x, u, 'o:', color = 'red', lw = 1.5, ms = 5, label = 'Initial')
#
####Periodic B.C's
#for n in range(nt):
#    un = u.copy()
#    for i in range(1, nx-1):
#        u[i] = un[i] - un[i]*(dt/dx)*(un[i] - un[i-1])\
#             + nu*(dt / dx**2)*(un[i+1] - 2*un[i] + un[i-1])
#        u[0] = un[0] - un[0]*(dt/dx)*(un[0] - un[-2])\
#             + nu*(dt / dx**2)*(un[1] - 2*un[0] + un[-2])
#        u[-1] = u[0] 
#u_ana = np.asarray([ufunc(nu, nt*dt, xi )for xi in x])
#
#plt.plot(x, u, 'o:', color = 'blue', lw = 1.5, ms = 5, label = 'Computational')
#plt.plot(x, u_ana, 'o:', color = 'purple', lw = 1.5, ms = 5, label = 'Analytical')
#plt.legend(loc = 'upper right', fontsize = 15)
#plt.show()

###Array operations with Numpy
import numpy as np
import timeit
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx)) ##create a 1xn vector of 1's
un = numpy.ones((ny, nx)) 

###Assign initial conditions

u[int(.5 / dy): int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

u = numpy.ones((ny, nx))
u[int(.5 / dy): int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

%%timeit
for n in range(nt + 1): #loop across number of time steps
    un = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c * dt / dx * 
                                  (un[j, i] - un[j, i - 1])) - 
                                  (c * dt / dy * 
                                   (un[j, i] - un[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1
















