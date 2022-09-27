# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 20:40:57 2022

@author: afdr9
"""
import numpy as np
import matplotlib.pyplot as plt

nx = 41
dx = 2/(nx-1)
nt = 25
dt = .025
c = 1

u = np.ones(nx)
u[int(.5/dx) : int(1/dx + 1)] = 2


plt.figure(figsize=(10,10))
plt.style.use(['bmh'])
plt.title('STEP 1', size = 15)
plt.xlabel('Distance', size = 15)
plt.ylabel('Velocity', size = 15)
plt.tick_params(axis = 'both', labelsize = 15)
plt.plot(np.linspace(0,2,nx), u, 'o:', color = 'blue', lw = 1.5, ms = 6, label = 'Initial value')


un = np.ones(nx)
for n in range(nt):
    un = u.copy()
    for i in range(1,nx):
        u[i] = un[i] - c * (dt/dx) * (un[i] - un[i-1])

plt.plot(np.linspace(0,2,nx), u, 'o:', color = 'red', lw = 1.5, ms = 6, label = 'Solution')
plt.legend(loc = 'upper right', fontsize = 15)
plt.show()
