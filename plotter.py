import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from particle_system import *

class Trajectory:
     def __init__(self):
          self.times = []
          self.positions = []

def _get_x_bounds(system: ParticleSystem): 
        positions = [p.position for p in system.particles]
        bound = max(abs(min(positions)),abs(max(positions)))
        bound *= 1.5
        return (-bound, bound)

def plot_ghost_state(system: ParticleSystem, totalTime: float, steps: int):
    fig, ax = plt.subplots()
    ax.set(
        xlim=(_get_x_bounds(system)),
        ylim=(0,totalTime), 
    )

    trajectories = compute_ghost_trajectories(system, totalTime, steps)
    for i, trajectory in enumerate(trajectories):
        ax.plot(trajectory.positions, 
                trajectory.times, 
                color=matplotlib.colors.hsv_to_rgb((i/len(trajectories),1,1)), 
                lw=0.5)
    plt.show()

def compute_ghost_trajectories(system: ParticleSystem, totalTime: float, steps: int) -> list[Trajectory]:
    trajectories = []
    for i, particle in enumerate(system.particles):
        new_trajectory = Trajectory()
        new_trajectory.times = np.linspace(0, totalTime, steps)
        new_trajectory.positions = [system.particles[i].evaluate_ghost_state(time) for time in new_trajectory.times]
        trajectories.append(new_trajectory)
    return trajectories
    
def plot_evolution(system: ParticleSystem, totalTime: float, steps: int):
    deltaTime = totalTime/steps
    particle_trajectories = compute_ghost_trajectories(system, totalTime)
         
