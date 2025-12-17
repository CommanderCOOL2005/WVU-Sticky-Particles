This python code simulates and graphs the evolution of a finite weighted sum of dirac masses under the Repulsive Pressureless Euler Poisson system.

# Background

The Repulsive Pressureless Euler Poisson System is defined to be

$\partial_t(\rho) + \partial_x(\rho v) = 0$ 

$\partial_t(\rho v) + \partial_x(\rho v^2) = \frac{1}{2}(\text{sgn} * \rho)\rho$

The code simulates and plots the paths of a weighted dirac sum, or if your initial rho distribution is equal to 

$\rho_0 = \sum_{i=1}^n m_i \delta_{y_i}$

# Usage
To use this code, import ```plotter.py``` into a python file. To create a particle system, the ```ParticleSystem()``` constructor must be called with an input of a list of particles created by the ```Particle()``` constructor. Note that passing in a list of particles into a ```ParticleSystem``` object will:
1. Divide the mass of each particle by the total mass of the list of particles so that the ```ParticleSystem```'s total mass will be one. 
2. Shift the position of each particle so that the center of mass of the system is zero.
3. Shift the velocity of each particle so that the total momentum of the system is zero.

Adjustments 2 and 3 can be done because the system is [Galilean Invariant](https://en.wikipedia.org/wiki/Galilean_invariance), and adjustment 1 is to make the measure a [probability measure](https://en.wikipedia.org/wiki/Probability_measure). This "normalization" can be overwritten if desired.

Additionally, the system will also set the acceleration of each particle to follow the Repulsive Pressureless Euler Poisson System. The exact acceleration set is $a_i=\frac{1}{2}\left(\sum_{i<j}m_j-\sum_{i>j}m_j\right)$, where $a_i$ denotes the acceleration of the i'th particle, and $m_j$ denotes the mass of the j'th particle. Note that the particles are ordered from left to right starting at i=0. This formula follows from the second part of RPEP. This can also be overwritten.

A perfect solution is one such that all collisions are "glancing" (i.e. at the precise moment of collision particles have the same velocity). A perfect solution is one where the velocity of each particle is continuous. If you want the velocities to be set such that the solution is perfect, you can call the ``assign_perfect_solution()`` method found in ``ParticleSystem`` once the system object has been created. 

Once a ```ParticleSystem``` has been created, you can simulate and graph the system using the ```plot_solution()``` method found in ```plotter.py```.

# Examples

Here are some examples:

**Example 1:** Plotting $p_1=(m_1,\frac{1}{2},y_1=-1,v_1=1),p_2(m_2=\frac{1}{2},y_1=1,v_1=-1)$
```python
particle_input = [Particle(0.5,-1,1),Particle(0.5,1,-1)]
system = ParticleSystem(particles=particle_input)
plot_solution(system, total_time=3, steps=200)
```
The resulting graph is:
![Graph of example 1](/example_figures/fig_1.png)


**Example 2:** Plotting 40 particles with random masses between 0.5 and 1, positions between -1 and 1, and velocities between -1 and 1
```python
particle_input = [Particle(uniform(0.5,1), uniform(-1,1), uniform(-1,1)) for i in range(40)]
system = ParticleSystem(particles=particle_input)
plot_solution(system, total_time=3, steps=200)
```

The resulting graph is:
![Graph of example 2](/example_figures/fig_2.png)

**Example 3:** Plotting 40 particles with random masses and positions and using ```assign_perfect_solution()``` to assign velocities such that the solution is perfect
```python
particle_input = [Particle(uniform(0.5,1), uniform(-1,1)) for i in range(40)]
system = ParticleSystem(particles=particle_input)
system.assign_perfect_solution() #assigns a perfect solution
plot_solution(system, total_time=3, steps=200)
```

The resulting graph is:
![Graph of example 3](/example_figures/fig_3.png)

**Example 4:** Same setup as example 3 except the ghost state of the system is plotted as well

```python
particle_input = [Particle(uniform(0.5,1), uniform(-1,1)) for i in range(40)]
system = ParticleSystem(particles=particle_input)
system.assign_perfect_solution() #assigns a perfect solution
plot_solution(system, total_time=3, steps=200, plot_ghost_state=True) #plot ghost state as well
```
The resulting graph is:
![Graph of example 4](/example_figures/fig_4.png)

# Dependencies
The main sticky particle solver uses matplotlib version 3.10.3. region_plotter.py used  pygame version 2.6.1.
