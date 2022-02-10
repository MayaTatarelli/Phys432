"""
Insert description of the code
@author : Insert your name
@collab: List the names of your collaborators if any Insert the date
"""

import numpy as np
import matplotlib . pyplot as pl

dt = 1 #insert appropriate timestep
Nsteps = 100 #insert appropriate total number of timesteps

## Setting up initial conditions (vortex centres and circulation) 
# Vortex rings
y_v = np.array([50.0,-50.0,50.0,-50.0])
x_v = np.array([-50.0,-50.0,50.0,50.0])
k_v = np.array ([-10.0,-10.0,-10.0,-10.0]) #insert the line vortex constant k of the 4 vortices 

# Setting up the plot
pl.ion()
fig , ax = pl.subplots(1,1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, 'k+', markersize=10)
#play around with the marker size and type as you see fit

# draw the initial velocity streamline
ngrid = 100 #insert the dimension of your simulation grid
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
#360j sets the resolution of the cartesian grid; play around with it as you see fit
vel_x = np.zeros(np.shape(Y)) #this holds x−velocity 
vel_y = np.zeros(np.shape(Y)) #this holds y−velocity

# masking radius for better visualization of the vortex centres 
r_mask = 2 #insert the radius of the mask around the vortex centres
#within this mask, you will not plot any streamline
#so that you can see more clearly the movement of the vortex centres
#R = np.sqrt(X**2+Y**2)

for i in range(len(x_v)): #looping over each vortex
    # insert lines for computing the total velocity field
    # insert lines for masking (set the masking area to NaN)
    vel_x += k_v[i]*(Y-y_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
    vel_y += -k_v[i]*(X-x_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
    #vel_x[R<r_mask] = np.nan
    #vel_y[R<r_mask] = np.nan

#/(2*np.pi)

# set up the boundaries of the simulation box 
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines 
ax.streamplot(X, Y, vel_x , vel_y , density=[1,1])
#play around with density as you see fit; 
#see the API documentation for more detail

fig.canvas.draw()

# Evolution
u_x = np.zeros(4)
u_y = np.zeros(4)
count = 0
while count < Nsteps:
    ## Compute and update advection velocity 
    # insert lines to re−initialize the total velocity field
    for i in range(len(x_v)):
        # insert lines to compute the total 
        # advection velocity on each vortex
        if i == 0:
            u_x[0] = -k_v[0]/(y_v[0] - y_v[1])
            u_y[0] = -k_v[0]/(x_v[0] - x_v[2]) 
            
        if i == 1:
            u_x[1] = -k_v[1]/(y_v[0] - y_v[1])
            u_y[1] = -k_v[1]/(x_v[3] - x_v[1])  
            
        if i == 2:
            u_x[2] = -k_v[2]/(y_v[2] - y_v[3])
            u_y[2] = -k_v[2]/(x_v[2] - x_v[0])  
            
        if i == 3:
            u_x[3] = -k_v[3]/(y_v[2] - y_v[3])
            u_y[3] = -k_v[3]/(x_v[1] - x_v[3])                        
        
    # insert lines to update the positions of vortices
    for i in range(len(x_v)):
        y_v[i] += u_y[i]*dt
        x_v[i] += u_x[i]*dt  
    
    #print(y_v)
    #print(x_v)
    #count += 1
    
    # insert lines to re−initialize the total velocity field 
    vel_x = vel_x*0.0
    vel_y = vel_y*0.0
    
    for i in range(len(x_v)):
        # insert lines to update the streamlines and masking
        vel_x += k_v[i]/(2*np.pi)*(Y-y_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
        vel_y += -k_v[i]/(2*np.pi)*(X-x_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
        
    ## update plot
        
    # the following two lines clear out the previous streamlines
    ax.collections = []
    ax.patches = []
    
    p.set_xdata(x_v) 
    p.set_ydata(y_v)
    
    ax.streamplot(X, Y, vel_x, vel_y, density=[1, 1])
    
    fig.canvas.draw()
    pl.pause(0.001) #play around with the delay time for better visualization
    count += 1
    

