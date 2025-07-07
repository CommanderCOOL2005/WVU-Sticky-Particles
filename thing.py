import copy

import matplotlib.pyplot as plt

from particle import *
from particle_system import *
from plotter import *

while True:
    initial_input = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(8)]
    initial_system = ParticleSystem(initial_input)
    initial_system.assign_perfect_solution()
    initial_system.print_particles()
    target_index = 0
    for i in (0.1,0.11,0.05,0.025):
        # initial_input = [Particle(1,(i - (n+1)/2)/((n-1)/2),-(i - (n+1)/2)/((n-1)/2)) for i in range(5)]
        step_size = i
        init_pos = initial_system.particles[target_index].position - 1
        final_pos = initial_system.particles[target_index].position + 1
        pos = init_pos
        init_vel = initial_system.particles[target_index].velocity - 1
        final_vel = initial_system.particles[target_index].velocity + 1
        vel = init_vel

        while vel <= final_vel:
            system_copy = copy.deepcopy(initial_system)
            system_copy.particles[target_index].velocity = vel
            pos = init_pos
            while pos <= final_pos:
                system_copy.particles[target_index].position = pos
                system_copy.normalize_system()
                is_eq = system_copy.is_equilibrium_solution()
                plt.plot(pos,vel,"s",color='green' if is_eq else 'red', ms=10)
                pos += step_size
            vel += step_size
        plt.show()