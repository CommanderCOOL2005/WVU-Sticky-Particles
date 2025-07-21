from random import *
from plotter import *

# An example main.py file

particle_input = [Particle(mass=uniform(0.5,1), position=uniform(-1,1)) for i in range(7)]
#Generate a random input of 7 particles with mass in the range of [0.5, 1] and position [-1, 1] 

system = ParticleSystem(particles=particle_input)
# Create a ParticleSystem object. Note that this object automatically adjusts the particles in particle input

system.assign_random_signed_velocities(0.5,1)
# Assign random velocities with magnitude from [0.5,1], which point right is particle is to the left of y=0, and point left is particle to the right of y=0 

plot_solution(
    system,                     # The system to plot
    plot_until_time_bound=True, # Plots until the maximum time of equilibrium
    steps=200,                  # Plots approximately 200 points
    plot_quadratic_bound=False, # Plots the quadratic bound (if true)
    plot_ghost_state=False      # Plots the ghost state (if true)
)
