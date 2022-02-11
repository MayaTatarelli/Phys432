"""
This is a simulation of the leapfrogging behaviour 
when two vortex rings interact with each other
@author : Maya Tatarelli
@collab: Yuliya Shpunarska
Feb 10th 2022
"""
#import statements
import numpy as np
import matplotlib . pyplot as pl

#with this timestep and number of steps, 
#it allows for a clear visual of the movement and interactions of the vortices
dt = 1 #insert appropriate timestep
Nsteps = 20 #insert appropriate total number of timesteps

## Setting up initial conditions (vortex centres and circulation) 
# Vortex rings
#The x-coord for the top and bottom of the left and right ring are the same, respectively
#index 0: top left, index 1: bottom left, index 2: top right, index 3: bottom right
y_v = np.array([35.0,-35.0,35.0,-35.0])
x_v = np.array([-20.0,-20.0,20.0,20.0])
k_v = np.array ([-10.0,10.0,-10.0,10.0]) #negative=ccw circulation, positive=cw circulation

# Setting up the plot
pl.ion()
fig , ax = pl.subplots(1,1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, 'k+', markersize=10)
#play around with the marker size and type as you see fit

# draw the initial velocity streamline
#This size allows for a clear unobstructed visual of the vortices
ngrid = 100 #insert the dimension of your simulation grid
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
#360j sets the resolution of the cartesian grid; play around with it as you see fit
vel_x = np.zeros(np.shape(Y)) #this holds x−velocity 
vel_y = np.zeros(np.shape(Y)) #this holds y−velocity

# masking radius for better visualization of the vortex centres 
r_mask = 2 #insert the radius of the mask around the vortex centres
#within this mask, you will not plot any streamline
#so that you can see more clearly the movement of the vortex centres

for i in range(len(x_v)): #looping over each vortex
    # insert lines for computing the total velocity field
    # insert lines for masking (set the masking area to NaN)
    #These compute the velocities at a given point in the grid due to each vortex
    #This comes from u = circulation/(2*pi*r) = k/r, where r is the distance between vortex and point
    #But it must be split into x and y components
    vel_x += k_v[i]*(Y-y_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2) #(k/r)*(delta y/r) to get x-comp of vel
    vel_y += -k_v[i]*(X-x_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2) #(k/r)*(delta x/r) to get y-comp of vel

# set up the boundaries of the simulation box 
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines 
ax.streamplot(X, Y, vel_x , vel_y , density=[1,1])
#play around with density as you see fit; 
#see the API documentation for more detail

fig.canvas.draw()

# Evolution
#These hold the x and y velocity components at a given vortex 
#due to the other three vortices
u_x = np.zeros(4)
u_y = np.zeros(4)
count = 0
while count < Nsteps:
    ## Compute and update advection velocity 
    # insert lines to re−initialize the total velocity field
    for i in range(len(x_v)):
        # insert lines to compute the total 
        # advection velocity on each vortex
        #similar to initialization of velocity field, 
        #these use the same formula to compute the advection velocity
        #but each vortex will skip itself in its calculation
        for j in range(len(x_v)):
            if(i != j):
                u_x[i] += k_v[j]*(y_v[i]-y_v[j])/((x_v[i]-x_v[j])**2+(y_v[i]-y_v[j])**2)
                u_y[i] += -k_v[j]*(x_v[i]-x_v[j])/((x_v[i]-x_v[j])**2+(y_v[i]-y_v[j])**2)
                  
        
    # insert lines to update the positions of vortices
    for i in range(len(x_v)):
        y_v[i] += u_y[i]*dt
        x_v[i] += u_x[i]*dt  
    
    # insert lines to re−initialize the total velocity field 
    vel_x = vel_x*0.0
    vel_y = vel_y*0.0
    
    for i in range(len(x_v)):
        # insert lines to update the streamlines and masking
        #This will update the streamlines in the same way they were initialized
        #But with the new vortex positions
        vel_x += k_v[i]*(Y-y_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
        vel_y += -k_v[i]*(X-x_v[i])/((X-x_v[i])**2+(Y-y_v[i])**2)
        
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
    

