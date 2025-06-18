import math

class Particles:
    def __init__(self, particles: List[Particle], rel_vel=0):
        self.particles = particles
        self.rel_vel = rel_vel  # relative velocity of the system

    def center(self):
        return sum(p.mass * p.position for p in self.particles) / sum(p.mass for p in self.particles)

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
