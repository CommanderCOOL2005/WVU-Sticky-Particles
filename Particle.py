import math 

class Particle:
    def __init__(self, mass, position):
        self.mass = mass
        self.position = position

    def __repr__(self):
        return f"Particle(m={self.mass}, y={self.position})"
    
    def __add__(self, other):
        if isinstance(other, Particle):
            return Particle(self.mass + other.mass, (self.mass * self.position + other.mass * other.position) / (self.mass + other.mass))
        raise TypeError("Can only add another Particle instance.")
    
    # time of perfect collision
    def perfect_time(self, other):
        if isinstance(other, Particle):
            return 2 * math.sqrt((other.position - self.position) / (self.mass + other.mass))
        raise TypeError("Can only compute time with another Particle instance.")
    
    # velocity difference, v1 - v2, needed for perfect solution.
    def perfect_difference(self, other):
        if isinstance(other, Particle):
            return math.sqrt((other.position - self.position) * (self.mass + other.mass))
        raise TypeError("Can only compute difference with another Particle instance.")
