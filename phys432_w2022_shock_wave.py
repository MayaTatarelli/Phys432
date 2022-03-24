"""
This is a simulation of an adibatic strong shock
using the donor cell advection scheme to solve the full
hydro equations

@author: Maya Tatarelli
March 24th 2022
"""
import numpy as np
import matplotlib.pyplot as pl

# Physical properties
adiabatic_index = 5.0/3.0

def get_pressure(f1,f2,f3):
    # computes the pressure for the Euler equations
    pressure = adiabatic_index * (f3 - 0.5 * f2 * f2 / f1)
    return pressure

# Set up the grid, time and grid spacing, and the sound speed squared
Ngrid = 250
Nsteps = 2500
dx = 100.0/np.float64(Ngrid)
courant_friedrichs_lewy_number = 0.05
dt = dx*courant_friedrichs_lewy_number/0.5 # CFL=c*dt/dx where c = 0.5 (advection speed)

u = np.zeros(Ngrid+1) # velocity
x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.zeros(Ngrid) # rho x e_tot

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update
    return f

# Apply initial Gaussian perturbation to total energy
Amp, sigma = 1.0, np.float64(Ngrid)/(10.0*(np.float64(Ngrid)/100.0))
f1 = 1.0*f1 # rho(t=0) = 1
f2 = 0.0*f2 # (rho x u) (t=0) = 0
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2) # (rho x e_tot) (t=0) = 1 + delta * (rho x e_tot)

# update mach number based on initialization
u_at_flux_points = f2/f1 # velocity
cs2 = adiabatic_index*get_pressure(f1,f2,f3)/f1 # sound speed square
mach_number = u_at_flux_points/np.sqrt(cs2) # mach number

# plotting
pl.ion()

fig, axs = pl.subplots(2)
fig.suptitle('Shock Simulation')
x1, = axs[0].plot(x, f1, 'ro')
x2, = axs[1].plot(x, mach_number, 'ro')

axs[0].set_xlim([0, dx*Ngrid+1])
axs[1].set_xlim([0, dx*Ngrid+1])
axs[0].set_ylim([0.0, 5.0])
axs[1].set_ylim([0.0, 1.2])

axs[0].set_xlabel('x')
axs[1].set_xlabel('x')
axs[0].set_ylabel('Density')
axs[1].set_ylabel('Mach Number')

fig.canvas.draw()

for ct in range(Nsteps):
    # (1) advection velocity at the cell interfaces
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))

    # (2) advect density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)

    # (3) add the source term
    # -- momentum equation
    f2[1:-1] = f2[1:-1] - 0.5 * (dt / dx) * (get_pressure(f1[2:],f2[2:],f3[2:]) - get_pressure(f1[:-2],f2[:-2],f3[:-2]))

    # (3) correct for source term at the boundary (reflective)
    # -- momentum equation
    f2[0] = f2[0] - 0.5 * (dt / dx) * (get_pressure(f1[1],f2[1],f3[1]) - get_pressure(f1[0],f2[0],f3[0]))
    f2[-1] = f2[-1] - 0.5 * (dt / dx) * (get_pressure(f1[-1],f2[-1],f3[-1]) - get_pressure(f1[-2],f2[-2],f3[-2]))

    # (4) re-calculate advection velocities with the inclusion of the boundaries
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))
    u[0] = 0.0 # for reflective boundary conditions
    u[-1] = 0.0 # for reflective boundary conditions

    # (5) advect energy
    f3 = advection(f3, u, dt, dx)

    # (6) add the source term
    # -- energy equation
    u_left = 0.5 * (u[2:-1] + u[3:])
    u_right = 0.5 * (u[:-3] + u[1:-2])
    f3[1:-1] = f3[1:-1] - 0.5 * (dt / dx) * (u_left * get_pressure(f1[2:],f2[2:],f3[2:]) - u_right * get_pressure(f1[:-2],f2[:-2],f3[:-2]))

    # (7) correct for source term at the boundary (reflective)
    # -- energy equation
    f3[0] = f3[0] - 0.5 * (dt / dx) * ( 0.5*(u[1]+u[2]) * get_pressure(f1[1],f2[1],f3[1]) - 0.5*(u[0]+u[1]) * get_pressure(f1[0],f2[0],f3[0]))
    f3[-1] = f3[-1] - 0.5 * (dt / dx) * ( 0.5*(u[-1]+u[-2]) * get_pressure(f1[-1],f2[-1],f3[-1]) - 0.5*(u[-2]+u[-3]) * get_pressure(f1[-2],f2[-2],f3[-2]))

    # (8) update sound speed
    cs2 = adiabatic_index*get_pressure(f1,f2,f3)/f1 # sound speed square
    mach_number = 0.5*(np.abs(u[:-1]+u[1:]))/np.sqrt(cs2) # mach number; absolute value to remove minus sign indicating direction

    # For checking post and pre shock values
    # if(ct==1239):
    #     print("Pre-shock density is %2.4f" % np.amax(f1))
    #     print("Pre-shock mach number is %2.4f" % np.amax(mach_number))
    #     print("Post-shock density is %2.4f" % f1[np.argmax(f1)+20])
    #     print("Post-shock mach number is %2.4f" % mach_number[np.argmax(mach_number)+20])
    #     print("Min density is %2.4f" % f1[np.argmin(f1)])
    #     print("Mach number associated with min density is %2.4f" % mach_number[np.argmin(f1)])

    # update the plot
    axs[0].set_title("Time Step: %i, dx = %2.3e, dt = %2.3e" % (ct,dx,dt))
    x1.set_ydata(f1)
    x2.set_ydata(mach_number)
    fig.canvas.draw()
    pl.pause(0.001)
