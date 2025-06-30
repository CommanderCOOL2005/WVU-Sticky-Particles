# Given an array of particles, (m_i, y_i),
# where m is the mass and the y is the position,
# find the perfect solution for the repulsive pressureless Euler-Poisson system.
# The perfect solution is defined as the solution where all collions are glancing
# (both particles have the same velocity) and an equllibrium is reached.
# Additionally, total momentum of the system is 0.
# Conservation of momentum holds

from random import *

from particle import *
from particle_system import *
from plotter import *

input = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(5)]

system = ParticleSystem(input)
system.assign_random_signed_velocities()
plot_evolution(system, 5, 100)

