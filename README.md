This python code simulates and graphs the evolution of a finite weighted sum of dirac masses under the Repulsive Pressureless Euler Poisson system.

# Usage
To use this code, import ```plotter.py``` into a python file. To create a particle system, the ```ParticleSystem()``` constructor must be called with an input of a list of particles created by the ```Particle()``` constructor. Note that passing in a list of particles into a ```ParticleSystem``` object will:
1. Divide the mass of each particle by the total mass of the list of particles so that the ```ParticleSystem```'s total mass will be one. 
2. Shift the position of each particle so that the center of mass of the system is zero.
3. Shift the velocity of each particle so that the total momentum of the system is zero.
These can all be done because the system is [Galilean Invariant](https://en.wikipedia.org/wiki/Galilean_invariance). 
If you want the velocities to be set such that the solution is perfect, you can call the ``assign_perfect_solution()`` method found in ``ParticleSystem`` once the system object has been created. 

Once a ```ParticleSystem``` has been created, you can simulate and graph the system using the ```plot_solution()``` method found in ```plotter.py```.

# Examples

Here are some examples:

Example 1: Plotting $p_1=(m_1,\frac{1}{2},y_1=-1,v_1=1),p_2(m_2=\frac{1}{2},y_1=1,v_1=-1)$
```python
particle_input = [Particle(0.5,-1,1),Particle(0.5,1,-1)]
system = ParticleSystem(particles=particle_input)
plot_solution(system, total_time=3, steps=200)
```
The resulting graph is:
![Graph of example 1](/example_figures/fig_1.png)

# Dependencies
This code uses matplotlib version 3.10.3