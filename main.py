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

input = [Particle(1,-1,1), Particle(1,1,-1), Particle(1,0,0)]

system = ParticleSystem(input)
system.do_everything()
system.print_particles()
system.step_fancy(3)
system.print_particles()

