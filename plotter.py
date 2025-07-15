import copy
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

def _compute_ghost_trajectories(system: ParticleSystem, total_time: float, steps: int, time_shift: float) -> list[Trajectory]:
    trajectories = []
    for i in range(len(system.particles)):
        new_trajectory = Trajectory(system.particles[i].color)
        new_trajectory.times = np.linspace(0, total_time, steps)
        new_trajectory.positions = [system.particles[i].evaluate_ghost_state(time) for time in new_trajectory.times]
        new_trajectory.times += np.full((steps,),time_shift)
        trajectories.append(new_trajectory)
    return trajectories

def plot_ghost_trajectories(system: ParticleSystem, total_time: float, steps: int, time_shift: float, _fig, _ax):
    trajectories = _compute_ghost_trajectories(system, total_time, steps, time_shift)
    for i, trajectory in enumerate(trajectories):
        _ax.plot(trajectory.positions, 
                trajectory.times, 
                color=trajectory.color, 
                lw=1)

def _plot_ghost_state(system: ParticleSystem, total_time: float, steps: int, fig, ax):
    plot_ghost_trajectories(system, total_time, steps, 0, fig, ax)
    
def _plot_real_state(system: ParticleSystem, total_time: float, steps: int, fig = None, ax = None):
    system_copy = copy.deepcopy(system)
    elapsed_time = 0
    temp_total_steps = 0
    while True:
        next_collision = system_copy.get_next_collision()
        delta_time = min(next_collision.time, total_time - elapsed_time)
        substeps = None
        substeps = max(2,int(steps*delta_time/(total_time-elapsed_time)))
        plot_ghost_trajectories(system_copy, delta_time, substeps, elapsed_time, fig, ax)
        system_copy.advance(delta_time, next_collision)
        if (total_time-elapsed_time) != 0:
            temp = int(steps*delta_time/(total_time-elapsed_time))
            steps -= temp
        elapsed_time += delta_time
        temp_total_steps += substeps
        
        if(elapsed_time >= total_time):
             break

def plot_solution(system: ParticleSystem, total_time: float, steps: int, plot_real_state: bool= True, plot_ghost_state: bool = False):
    """Generates a plot of the system over time.

    Args:
        system (ParticleSystem): The system to plot
        total_time (float): The total time to plot the system for
        steps (int): Approximately the amount of points to plot for each particle (this will slightly vary due to how the system plots, but it will be very close to this number) 
        plot_real_state (bool, optional): Whether to plot the actual state of the system. Defaults to True.
        plot_ghost_state (bool, optional): Whether to plot the "ghost state" of the system, which is how the system would evolve if particles did not stick together and instead just passed through each other. Defaults to False.
    """
    axNum = int(plot_ghost_state) + int(plot_real_state)
    fig, axes = plt.subplots(axNum, figsize=(10,10))
    if axNum == 1: #this is the dumbest thing i have ever had to write why is matplotlib the stupidest library on earth why would it do something like this
        axes = [axes]
    xbounds = _get_x_bounds(system)
    for ax in axes:
        ax.set(xlim = xbounds,
               ybound = (0, total_time)) 
    axIndex = 0
    if plot_real_state:
        axes[axIndex].set(title = "Real State")
        _plot_real_state(system, total_time, steps, fig, axes[axIndex])
        axIndex += 1
    if plot_ghost_state:
        axes[axIndex].set(title = "Ghost State")
        _plot_ghost_state(system, total_time, steps, fig, axes[axIndex])
        axIndex += 1
    plt.show()

         
