from random import *

from plotter import *

particle_input = [Particle(0.5,-1,1),Particle(0.5,1,-1)]
system = ParticleSystem(particles=particle_input)
plot_solution(system, total_time=3, steps=200)
