import copy
from math import *
from random import *
from particle import Particle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from helper import *

class Collision:
    def __init__(self, time, indices = [], sizes = []):
        self.time = time
        self.indices = indices
        self.sizes = sizes
        pass

class ParticleSystem:
    def __init__(
                self, 
                particles: list[Particle], 
                override_normalization = False,
                override_acceleration = False,
                relative_velocity=0
                ):
        if particles is None:
            self.particles = []
        else:
            self.particles = sorted(particles, key=lambda particle: particle.position)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
        self.make_rgb_colors()
        self.rel_vel = relative_velocity  # relative velocity of the system

    def get_total_mass(self):
        return sum(p.mass for p in self.particles)
    
    def get_center_of_mass(self):
        return sum(p.mass * p.position for p in self.particles) / self.get_total_mass()
        
    def get_total_momentum(self):
        return sum(p.mass * p.velocity for p in self.particles)

    def set_accelerations(self):
        for i in range(len(self.particles)):
            new_accel = 0
            for j in range(len(self.particles)):
                if i == j:
                    continue
                elif i < j:
                    new_accel -= self.particles[j].mass
                else:
                    new_accel += self.particles[j].mass
            self.particles[i].acceleration = 0.5*new_accel

    def normalize_mass(self):
        total_mass = self.get_total_mass()
        for p in self.particles:
            p.mass /= total_mass

    def center_mass(self):
        com = self.get_center_of_mass()
        for p in self.particles:
            p.position -= com

    def center_momentum(self):
        total_momentum = sum(p.mass * p.velocity for p in self.particles)
        total_mass = sum(p.mass for p in self.particles)
        mean_velocity = total_momentum / total_mass

        for p in self.particles:
            p.velocity -= mean_velocity
    
    def configure(self, normalize_mass = False, center_mass = False, center_momentum = False):
        if normalize_mass:
            self.normalize_mass()
        if center_mass:
            self.center_mass()
        if center_momentum:
            self.center_momentum()

    def normalize_system(self, normalize_mass = True, center_mass = True, center_momentum = True):
        self.configure(normalize_mass, center_mass, center_momentum)
    
    def add_particle(self, particle: Particle, override_normalization = False, override_acceleration = False):
        self.particles.append(particle)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
    
    def add_particles(self, particles: list[Particle], override_normalization = False, override_acceleration = False):
        self.particles.extend(particles)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
    
    def get_next_collision(self):
        earliest = Collision(inf)
        lastParticleHadCollision = False
        for i in range(len(self.particles) - 1):
            p1 = self.particles[i]
            p2 = self.particles[i+1]
            discriminant = 2*p1.acceleration*(p2.position-p1.position) + 2*p2.acceleration*(p1.position - p2.position) + (p1.velocity-p2.velocity)**2
            leeway = -1e-15 # prevents very slight numerical error
            if leeway < discriminant < 0:
                discriminant = 0
            if discriminant < 0: #no collisions whatsoever, skip
                lastParticleHadCollision = False
                continue
            else:
                sqrtDisc = sqrt(discriminant)
                greaterTime = (p2.velocity - p1.velocity + sqrtDisc)/(p1.acceleration - p2.acceleration)
                lesserTime = (p2.velocity - p1.velocity - sqrtDisc)/(p1.acceleration - p2.acceleration)
                collisionTime = -13490091340
                if lesserTime > greaterTime:
                    temp = greaterTime
                    greaterTime = lesserTime
                    lesserTime = temp
                if greaterTime < 0: # both collision times are negative, skip
                    lastParticleHadCollision = False
                    continue
                if lesserTime < 0:  # one collision is negative, pick the other time
                    collisionTime = greaterTime
                else:               # pick the smallest time as both are positive
                    collisionTime = lesserTime
                
                if collisionTime > earliest.time:
                    lastParticleHadCollision = False
                    continue
                
                if collisionTime == earliest.time:        # add it to an existing collision
                    if lastParticleHadCollision:
                        earliest.sizes[-1] += 1
                    else:
                        earliest.indices.append(i)
                        earliest.sizes.append(2)
                        lastParticleHadCollision = True
                elif 0 < collisionTime < earliest.time:       # found new collision candidate, exclude zero
                    earliest.time = collisionTime
                    earliest.indices.clear()
                    earliest.sizes.clear()
                    earliest.indices.append(i)
                    earliest.sizes.append(2)
                    lastParticleHadCollision = True
                else:
                    continue
        return earliest

    def advance(self, totalTime: float, _collision : Collision = None):
        if totalTime <= 0:
            return
        if _collision == None:
            collision= self.get_next_collision()
        else:
            collision = _collision
        
        steppingTime = min(totalTime, collision.time)
        for particle in self.particles:
            particle.step_ghost_state(steppingTime)

        for i, collisionIndex in enumerate(collision.indices): #only runs if a collision will happen
            indexStart = collisionIndex
            indexEnd = collisionIndex+collision.sizes[i]

            particlesInCollision = self.particles[indexStart:indexEnd]
            newMass = sum(p.mass for p in particlesInCollision)
            newPos = particlesInCollision[0].position
            newVel = sum(p.mass*p.velocity for p in particlesInCollision)/newMass
            biggestMass = 0
            biggestParticle = None
            for p in particlesInCollision:
                if p.mass > biggestMass:
                    biggestParticle = p
            newColor = biggestParticle.color

            newParticle = Particle(newMass, newPos, newVel, color=newColor)

            del self.particles[indexStart:indexEnd]
            self.particles.insert(indexStart, newParticle)
            self.set_accelerations()

            self.advance(totalTime-collision.time)
            return
        return

    def advance_to_next_collision(self):
        collision = self.get_next_collision()
        time = 0
        if collision.time != inf:
            time = collision.time
        self.advance(time, collision)
        return time
    
    def is_equilibrium_solution(self):
        system_copy = copy.deepcopy(self)
        while True:
            time = system_copy.advance_to_next_collision()
            if time == 0:
                if len(system_copy.particles) > 1:
                    return False
                else:
                    return True
        

    # Find the time of the next perfect collision in the system.
    # Each collision gives us the following information:
    # v1 - v2 = sqrt((y2 - y1) * (m1 + m2))
    # After the collision, the particles will combined into one larger particle.
    def next_collision(self):
        min_time = float('inf')
        ghost = None
        k = None
        for i in range(len(self.particles) - 1):
            time = self.particles[i].perfect_time(self.particles[i+1])
            if time < min_time:
                min_time = time
                ghost = self.particles[i] + self.particles[i+1]
                k = i
        return min_time, ghost, k
    
    def get_perfect_solution(self):
        if len(self.particles) < 1:
            return []
        elif len(self.particles) == 1:
            return [0]
        elif len(self.particles) == 2:
            p1, p2 = self.particles
            return [p2.mass * sqrt((p2.position - p1.position)/(p1.mass + p2.mass)), (-1) * (p1.mass) * sqrt((p2.position - p1.position)/(p1.mass + p2.mass))]
        else:
            min_time, ghost, k = self.next_collision()
            p1, p2 = self.particles[k], self.particles[k+1]
            dif = p1.perfect_difference(p2)
            merged_system = ParticleSystem(self.particles[:k] + [ghost] + self.particles[k+2:])
            merged_solution = merged_system.get_perfect_solution()
            return merged_solution[:k] + [merged_solution[k] + (p2.mass * dif)/(ghost.mass), merged_solution[k] - (p1.mass * dif)/(ghost.mass)] + merged_solution[k+1:]

    def assign_perfect_solution(self):
        velocities = self.get_perfect_solution()
        for i in range(len(self.particles)):
            self.particles[i].velocity = velocities[i]
        self.center_momentum()

    def assign_total_question_mark_solution(self):
        for p in self.particles:
            p.velocity = -sgn(p.position)*sqrt(abs(2*p.position*p.acceleration))
        self.center_momentum()

    def assign_random_signed_velocities(self, a: float = 0, b:float = 1):
        for i in range(len(self.particles)):
            self.particles[i].velocity = (1 if self.particles[i].position < 0 else -1)*uniform(a, b)
        self.center_momentum()

    def make_rgb_colors(self):
        for i, particle in enumerate(self.particles):
            particle.color = matplotlib.colors.hsv_to_rgb((i/len(self.particles),1,0.75))

    def print_info(self):
        print(f"CHARACTERISTICS:\n" +
              f"{"Total Mass:":27} {round(self.get_total_mass(),5)}\n" +
              f"{"Center of Mass:":27} {round(self.get_center_of_mass(),5)}\n" +
              f"{"Total Momentum:":27} {round(self.get_total_momentum(),5)}\n\n" +
              f"DEBUG:\n" +
              f"{"Accelerations Set:":27} {self.flags["accelerations_set"]}\n" +
              f"{"Normalized Mass:":27} {self.flags["has_normal_mass"]}\n" +
              f"{"Zero Center of Mass:":27} {self.flags["has_zero_com"]}\n" +
              f"{"Zero Total Momentum:":27} {self.flags["has_zero_mom"]}\n"
            )
    
    def print_particles(self):
        precision = 5
        print(f"{"Mass:":13}{[round(p.mass, precision) for p in self.particles]}\n" +
              f"{"Position:":13}{[round(p.position, precision) for p in self.particles]}\n"+
              f"{"Velocity:":13}{[round(p.velocity, precision) for p in self.particles]}\n"+
              f"{"Acceleration:":13}{[round(p.acceleration, precision) for p in self.particles]}")
    
    def dump(self):
        for i in range(len(self.particles)):
            print(self.particles[i], end=',')
        print(self.particles[-1])