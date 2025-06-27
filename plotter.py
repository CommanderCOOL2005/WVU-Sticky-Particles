import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from particle_system import *

def get_x_bounds(system: ParticleSystem): 
        positions = [p.position for p in system.particles]
        bound = max(abs(min(positions)),abs(max(positions)))
        bound *= 1.5
        return (-bound, bound)

def plot_ghost_state(system: ParticleSystem):
    height = 15
    fig, ax = plt.subplots()
    ax.set(
        xlim=(get_x_bounds(system)),
        ylim=(0,height), 
        xlabel=r'Position'
    )

    for i in range(len(system.particles)):
        y = np.linspace(0,height,100)
        x = [system.particles[i].evaluate_ghost_state(time) for time in y]
        ax.plot(x,y, color=matplotlib.colors.hsv_to_rgb((i/len(system.particles),1,1)), lw=0.5)
    plt.show()

def plot_real_state(system: ParticleSystem, end_time: float, steps: int):
    deltaTime = end_time/steps
    fig, ax = plt.subplots()
    ax.set(
        xlim=(get_x_bounds(system)),
        ylim=(0,end_time+1)
    )
    particle_histories = {} 
    for p in system.particles:
        particle_histories[p.id] = [(0, p.position)]
    for i in range(steps):
        step_real(system, deltaTime)
        current_time = i * deltaTime
        for p in system.particles:
            if p.id not in particle_histories:
                particle_histories[p.id] = []
            particle_histories[p.id].append((current_time,p.position))
    for history in particle_histories.values():
        times, positions = zip(*history)
        ax.plot(positions, times, lw=0.5)
    plt.show()

