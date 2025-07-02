from particle import *
from particle_system import *
from plotter import *

n = 5
initial_input = [Particle(1,(i - (n+1)/2)/((n-1)/2),-(i - (n+1)/2)/((n-1)/2)) for i in range(5)]
target_index = 0
system = ParticleSystem(initial_input)
plot(system, 3, 200)