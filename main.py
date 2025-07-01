# Given an array of particles, (m_i, y_i),
# where m is the mass and the y is the position,
# find the perfect solution for the repulsive pressureless Euler-Poisson system.
# The perfect solution is defined as the solution where all collions are glancing
# (both particles have the same velocity) and an equllibrium is reached.
# Additionally, total momentum of the system is 0.
# Conservation of momentum holds

# TODO: method to say if given configuration converges
# TODO: Brute force to see if certain configurations are possible (ie in 5 particles the first two collide then the next two then they all collide together)
# TODO: Make UI better
# TODO: Make graphs of other things (N_t, V_t, Energy)
# TODO: Mark collisions

from random import *

from particle import *
from particle_system import *
from plotter import *

innputte = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(4)]
system = ParticleSystem(innputte)
system.assign_total_question_mark_solution()
plot(system, 3, 200)
