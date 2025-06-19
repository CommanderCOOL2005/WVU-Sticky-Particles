# Given an array of particles, (m_i, y_i),
# where m is the mass and the y is the position,
# find the perfect solution for the repulsive pressureless Euler-Poisson system.
# The perfect solution is defined as the solution where all collions are glancing
# (both particles have the same velocity) and an equllibrium is reached.
# Additionally, total momentum of the system is 0.
# Conservation of momentum holds.

# Goal: return array of velocities v_i such that the system becomes a perfect solution.
# User inputs the array of tuples in ascending order of y.

import math
from typing import List, Tuple
from particle import *
from particle_system import *
import random as r
from plotter import *

# system = ParticleSystem()
# for i in range(20):
#     particle = Particle(r.uniform(0.5,1), r.uniform(-1,1), r.uniform(-1,1))
#     system.add_particle(particle)
# system.calibrate()
# plot_ghost_state(system)

# particles_input = input("Enter particles as (m_i, y_i) tuples in ascending order of y: ")
# # Ensure the input is a valid list of tuples
# try:
#     particles_input = eval(particles_input)
# except (SyntaxError, NameError):
#     raise ValueError("Invalid input format. Please enter a list of tuples in the form [(m1, y1), (m2, y2), ...].")
# # Validate that the input is a list of tuples
# if not isinstance(particles_input, list) or not all(isinstance(p, tuple) and len(p) == 2 for p in particles_input):
#     raise ValueError("Input must be a list of tuples in the form [(m1, y1), (m2, y2), ...].")
# # Validate that the masses are positive and positions are numeric
# for m, y in particles_input:
#     if not isinstance(m, (int, float)) or not isinstance(y, (int, float)):
#         raise ValueError("Each particle must be a tuple of (mass, position) with numeric values.")
#     if m <= 0:
#         raise ValueError("Mass must be positive.")
# # Ensure that the particles are in ascending order of y
# def is_ascending_order(particles):
#     return all(particles[i][1] < particles[i + 1][1] for i in range(len(particles) - 1))
# if not is_ascending_order(particles_input):
#     raise ValueError("Particles must be in ascending order of y.")

# input = [Particle(m, y) for m, y in particles_input]

input = [
    Particle(1, -2),
    Particle(2, -1),
    Particle(3, 1.5),
    Particle(4, 2.0)
]

system = ParticleSystem(input)
print("Solution:" + str(system.perfect_solution()))

# Example: Create a system with particles that have initial velocities
print("\nDemonstrating particle system evolution with collisions:")

# Create particles with positions and velocities
particles_with_velocities = [
    Particle(mass=1.0, position=-3.0, velocity=2.0),   # Moving right
    Particle(mass=1.5, position=-1.0, velocity=1.0),   # Moving right slower
    Particle(mass=2.0, position=1.0, velocity=-1.5),   # Moving left
    Particle(mass=1.2, position=3.0, velocity=-0.8),   # Moving left slower
]

# Create system and plot evolution
system_with_velocities = ParticleSystem(particles_with_velocities)
plot_system_evolution(system_with_velocities, max_time=5.0)