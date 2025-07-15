import math 
import matplotlib

class Particle:
    _next_id = 0
    def __init__(self, mass, position, velocity=0, acceleration=0, color=(1,0,0)):
        """A particle representing one weighted dirac mass in the Repulsive Pressureless Euler Poisson system.

        Args:
            mass (_type_): The mass of the particle
            position (_type_): The position of the particle
            velocity (int, optional): The velocity of the particle. Note that this value *might* be overwritten by ``ParticleSystem`` if certain methods within ``ParticleSystem`` are called. Defaults to 0.
            acceleration (int, optional): The acceleration of the particle. Note that this value *will* be overwritten when passed into a ``ParticleSystem`` object unless you override that behavior within the ``ParticleSystem`` object. Defaults to 0.
            color (tuple, optional): The color this particle will have when plotted. Defaults to (1,0,0).
        """
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.id = Particle._next_id
        self.color = color
        Particle._next_id += 1

    def __repr__(self):
        return f"Particle(mass={self.mass}, position={self.position}, velocity={self.velocity}, acceleration={self.acceleration})"
    
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
    
    def step_ghost_state(self, time):
        self.position += self.velocity * time + 0.5 * self.acceleration * time * time
        self.velocity += self.acceleration * time
