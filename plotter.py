import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from particle_system import *

class Trajectory:
     def __init__(self, color):
          self.times = []
          self.positions = []
          self.color = color

def _get_x_bounds(system: ParticleSystem): 
        positions = [p.position for p in system.particles]
        bound = max(abs(min(positions)),abs(max(positions)))
        bound *= 1.5
        return (-bound, bound)

def compute_ghost_trajectories(system: ParticleSystem, total_time: float, steps: int, time_shift: float = 0) -> list[Trajectory]:
    trajectories = []
    for i in range(len(system.particles)):
        new_trajectory = Trajectory(system.particles[i].color)
        new_trajectory.times = np.linspace(0, total_time, steps)
        new_trajectory.positions = [system.particles[i].evaluate_ghost_state(time) for time in new_trajectory.times]
        new_trajectory.times += np.full((steps,),time_shift)
        trajectories.append(new_trajectory)
    return trajectories

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
                color=trajectory.color, 
                lw=1)
    # plt.show()
    
def plot_evolution(system: ParticleSystem, total_time: float, steps: int):
    elapsed_time = 0
    fig, ax = plt.subplots()
    ax.set(
        xlim=(_get_x_bounds(system)),
        ylim=(0,total_time), 
    )
    while True:
        next_collision = system.get_next_collision()
        delta_time = min(next_collision.time, total_time - elapsed_time)
        substeps = max(2,int(steps*delta_time/(total_time-elapsed_time)))
        trajectories = compute_ghost_trajectories(system, delta_time, substeps, elapsed_time)
        for i, trajectory in enumerate(trajectories):
            ax.plot(trajectory.positions, 
                    trajectory.times,
                    color=trajectory.color,
                    lw=1) 
            #color=matplotlib.colors.hsv_to_rgb((i/len(trajectories),1,1))
        system.advance(delta_time, next_collision)
        elapsed_time += delta_time
        steps -= substeps
        if(elapsed_time >= total_time):
             break
    plt.show()
    

         
