# Given an array of particles, (m_i, y_i),
# where m is the mass and the y is the position,
# find the perfect solution for the repulsive pressureless Euler-Poisson system.
# The perfect solution is defined as the solution where all collions are glancing
# (both particles have the same velocity) and an equllibrium is reached.
# Additionally, total momentum of the system is 0.
# Conservation of momentum holds.

# Goal: return array of velocities v_i such that the system becomes a perfect solution.
# User inputs the array of tuples in ascending order of y.

import math
from typing import List, Tuple
from Particle import *
from Particles import *

particles_input = input("Enter particles as (m_i, y_i) tuples in ascending order of y: ")
# Ensure the input is a valid list of tuples
try:
    particles_input = eval(particles_input)
except (SyntaxError, NameError):
    raise ValueError("Invalid input format. Please enter a list of tuples in the form [(m1, y1), (m2, y2), ...].")
# Validate that the input is a list of tuples
if not isinstance(particles_input, list) or not all(isinstance(p, tuple) and len(p) == 2 for p in particles_input):
    raise ValueError("Input must be a list of tuples in the form [(m1, y1), (m2, y2), ...].")
# Validate that the masses are positive and positions are numeric
for m, y in particles_input:
    if not isinstance(m, (int, float)) or not isinstance(y, (int, float)):
        raise ValueError("Each particle must be a tuple of (mass, position) with numeric values.")
    if m <= 0:
        raise ValueError("Mass must be positive.")
# Ensure that the particles are in ascending order of y
def is_ascending_order(particles):
    return all(particles[i][1] < particles[i + 1][1] for i in range(len(particles) - 1))
if not is_ascending_order(particles_input):
    raise ValueError("Particles must be in ascending order of y.")

#interpret particles as array of Particle objects
def interpret_particles(particles: List[Tuple[float, float]]) -> List['Particle']:
    return [Particle(m, y) for m, y in particles]