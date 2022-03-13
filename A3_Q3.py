"""
This is a simulation of lava flowing down an inclined slope 
@author : Maya Tatarelli
@collab: Yuliya Shpunarska
March 12th 2022
"""
#import statements
import numpy as np
import matplotlib . pyplot as pl

# Set up the grid, timestep and number of steps and advection and diffusion parameters
#with this timestep and number of steps, 
#it allows for a clear visual of the time-evolution of the velocity field of lava flow
Ngrid = 50
Nsteps = 1250
dt = 0.001
dx = 1

#Constants
#Units in cm, s
g = 1000 #cm/s^2 gravitational acceleration
D1 = 2.5E3 #cm^2/s viscosity (diffusion coefficient)
theta = np.pi/4 #radians
rho = 3E3 #density
beta1 = D1*dt/dx**2

x = np.arange(0, Ngrid*1., dx) / Ngrid # multiplying by 1. to make sure this is an array of floats not integers

f1 = np.copy(x)

# Set up plot
pl.ion()
fig, axes = pl.subplots(1,1)
axes.set_title('Simulation of lava flowing down an inclined slope')
axes.set_xlabel('Flow velocity')
axes.set_ylabel('Height from inclined slope')

# Plotting steady state in the background for reference
vel = (-g/D1)*np.sin(theta)*(0.5*x**2 - 1.0*x)
axes.plot(x, vel, 'b-')

# We will be updating these plotting objects
plt1, = axes.plot(x, f1, 'ro')

# Setting the axis limits for visualization
axes.set_xlim([0,1])
axes.set_ylim([0,0.7])

# this draws the objects on the plot
fig.canvas.draw()

#Time evolution of lava flow
for count in range(Nsteps):
    
    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A1 = np.eye(Ngrid) * (1.0 + 2.0 * beta1) + np.eye(Ngrid, k=1) * -beta1 + np.eye(Ngrid, k=-1) * -beta1

    ## Boundary conditions to keep the first element fixed
    # This ensures f in the first cell stays fixed at all times under diffusion
    A1[0][0] = 1.0

    A1[0][1] = 0
    
    # Stress-free boundary condition on the right
    A1[Ngrid - 1][Ngrid - 1] = 1.0 + beta1
    
    #Adding constant to account for gravity (from momentum equation)
    f1 += dt*g*np.sin(theta)/rho
    f1[0] = 0
    
    # Solving for the next timestep
    f1 = np.linalg.solve(A1, f1)
    
    # update the plot
    plt1.set_ydata(f1)
    
    #Draw new points
    fig.canvas.draw()
    pl.pause(0.001)