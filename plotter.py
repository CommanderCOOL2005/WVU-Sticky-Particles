import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from particle_system import *

def get_x_bounds(system: ParticleSystem): 
        candidates = [abs(-p.velocity/p.acceleration) for p in system.particles]
        greatest_y_bound = max(candidates)
        index = candidates.index(greatest_y_bound)
        bound = system.particles[index].evaluate_ghost_state(greatest_y_bound)
        return (-bound, bound)

def plot_ghost_state(system: ParticleSystem):
    system.check_calibration()
    fig, ax = plt.subplots()
    ax.set(
        xlim=(system.get_x_bounds()),
        ylim=(0,5), 
        xlabel=r'Position'
    )

    for i in range(len(system.particles)):
        y = np.linspace(0,10,100)
        x = [system.particles[i].evaluate_ghost_state(time) for time in y]

        ax.plot(x,y, color=matplotlib.colors.hsv_to_rgb((i/len(system.particles),1,1)), lw=0.5)
    plt.show()