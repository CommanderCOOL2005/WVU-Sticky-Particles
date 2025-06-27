import math 

class Particle:
    _next_id = 0
    def __init__(self, mass, position, velocity=0, acceleration=0):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.id = Particle._next_id
        Particle._next_id += 1

    def __repr__(self):
        return f"Particle(m={self.mass}, y={self.position}, v={self.velocity}, a={self.acceleration})"
    
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
    
    def evaluate_ghost_state(self, time):
        return self.position + self.velocity * time + 0.5 * self.acceleration * time * time
    
    def step(self, deltaTime):
        self.velocity += self.acceleration * deltaTime
        self.position += self.velocity * deltaTime
