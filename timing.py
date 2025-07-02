from random import *
import time

from particle import *
from particle_system import *
from plotter import *

average_time = 0
last_time = 0
number_of_trials = 100
for n in range(1,101):
    for i in range(number_of_trials):
        time_start = time.time()
        innputte = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(n)]
        system = ParticleSystem(innputte)
        system.assign_random_signed_velocities(1,2)
        is_eq = system.is_equilibrium_solution()
        time_end = time.time()
        average_time += time_end - time_start
    average_time /= number_of_trials
    print(f"n={n:5} | {round(1000*(average_time),2):10} | {round(1000*(average_time-last_time),2):10}")
    last_time = average_time
    # plot(system, 3, 200)