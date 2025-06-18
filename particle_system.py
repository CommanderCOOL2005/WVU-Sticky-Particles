import math
from particle import Particle

class ParticleSystem:
    def __init__(self, particles: list[Particle] = [], relative_velocity=0):
        self.particles = particles
        self.calibrated = False
        self.rel_vel = relative_velocity  # relative velocity of the system

    def get_total_mass(self):
        return sum(p.mass for p in self.particles)
    
    def get_center_of_mass(self):
        return (sum(particle.mass * particle.position for particle in self.particles) / self.get_total_mass(self))
    
    def add_particle(self, particle: Particle):
        self.particles.append(particle)
        self.calibrated = False
    
    def calibrate(self):
        total_mass = sum(particle.mass for particle in self.particles)
        for particle in self.particles:
            particle.mass /= total_mass
        for i in range(len(self.particles)):
            acceleration = 0
            for j in range(len(self.particles)):
                if i == j:
                    continue
                if i < j:
                    acceleration -= self.particles[j].mass
                else:
                    acceleration += self.particles[j].mass
            self.particles[i].acceleration = acceleration
        self.calibrated = True

    def check_calibration(self):
        if not self.calibrated:
            print("WARNING: Particle system is not calibrated. Calibrate if this was not your intention")

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
            return [math.sqrt((p2.position - p1.position)(p1.mass + p2.mass)) + self.rel_vel, self.rel_vel]
