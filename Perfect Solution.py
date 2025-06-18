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


class Particle:
    def __init__(self, m, y):
        self.m = m
        self.y = y

    def __repr__(self):
        return f"Particle(m={self.m}, y={self.y})"
    
    def __add__(self, other):
        if isinstance(other, Particle):
            return Particle(self.m + other.m, (self.m * self.y + other.m * other.y) / (self.m + other.m))
        raise TypeError("Can only add another Particle instance.")
    
    # time of perfect collision
    def perfect_time(self, other):
        if isinstance(other, Particle):
            return 2 * math.sqrt((other.y - self.y) / (self.m + other.m))
        raise TypeError("Can only compute time with another Particle instance.")
    
    # velocity difference, v1 - v2, needed for perfect solution.
    def perfect_difference(self, other):
        if isinstance(other, Particle):
            return math.sqrt((other.y - self.y) * (self.m + other.m))
        raise TypeError("Can only compute difference with another Particle instance.")

class Particles:
    def __init__(self, particles: List[Particle], rel_vel=0):
        self.particles = particles
        self.rel_vel = rel_vel  # relative velocity of the system

    def center(self):
        return sum(p.m * p.y for p in self.particles) / sum(p.m for p in self.particles)

    # Find the time of the next perfect collision in the system.
    # Each collision gives us the following information:
    # v1 - v2 = sqrt((y2 - y1) * (m1 + m2))
    # After the collision, the particles will combined into one larger particle.
    def next_collision(self):
        min_time = float('inf')
        pair = None
        for i in range(len(self.particles) - 1):
            time = self.particles[i].perfect_time(self.particles[i+1])
            if time < min_time:
                min_time = time
                pair = (self.particles[i], self.particles[i+1])
                k = i
        return min_time, pair, k
    
    def perfect_solution(self):
        if len(self.particles) < 1:
            return []
        elif len(self.particles) == 1:
            return [0]
        elif len(self.particles) == 2:
            p1, p2 = self.particles
            return [math.sqrt((p2.y - p1.y)(p1.m + p2.m)) + self.rel_vel, self.rel_vel]
        else:
