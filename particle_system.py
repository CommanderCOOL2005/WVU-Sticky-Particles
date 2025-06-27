from math import *
from particle import Particle
import matplotlib.pyplot as plt
import numpy as np



class ParticleSystem:
    def __init__(self, particles: list[Particle] = None, relative_velocity=0):
        if particles is None:
            self.particles = []
        else:
            self.particles = sorted(particles, key=lambda particle: particle.position)
        self.total_mass = None
        self.center_of_mass = None
        self.total_momentum = None
        self.flags = {
            "total_mass_set": False,
            "center_of_mass_set": False,
            "total_momentum_set": False,
            "accelerations_set" : False,
            "has_normal_mass" : False,
            "has_zero_com" : False, #com = center of mass
            "has_zero_mom" : False  #mom = momentum
        }
        self.set_basics()
        self.rel_vel = relative_velocity  # relative velocity of the system

    def get_total_mass(self):
        if not self.flags["total_mass_set"]:
            self.total_mass = sum(p.mass for p in self.particles)
            self.flags["total_mass"] = True
        return self.total_mass
    
    def get_center_of_mass(self):
        if not self.flags["center_of_mass_set"]:
            self.center_of_mass = sum(p.mass * p.position for p in self.particles) / self.get_total_mass()
            self.flags["center_of_mass"] = True
        return self.center_of_mass
        
    def get_total_momentum(self):
        if not self.flags["total_momentum_set"]:
            self.total_momentum = sum(p.mass * p.velocity for p in self.particles)
            self.flags["total_momentum"] = True
        return self.total_momentum

    def set_accelerations(self):
        if not self.flags["accelerations_set"]:
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
            self.flags["accelerations_set"] = True

    def normalize_mass(self):
        if not self.flags["has_normal_mass"]:
            total_mass = self.get_total_mass()
            for p in self.particles:
                p.mass /= total_mass
            self.flags["has_normal_mass"] = True

    def center_mass(self):
        if not self.flags["has_zero_com"]:
            com = self.get_center_of_mass()
            for p in self.particles:
                p.position -= com
            self.flags["has_zero_com"] = True

    def center_momentum(self):
        if not self.flags["has_zero_mom"]:
            total_velocity = sum([p.velocity for p in self.particles])
            for p in self.particles:
                p.velocity -= total_velocity
            self.flags["has_zero_mom"] = True

    def set_basics(self):
        self.get_center_of_mass()
        self.get_total_mass()
        self.get_total_momentum()
    
    def do_things(self, set_accelerations = False, normalize_mass = False, center_mass = False, center_momentum = False):
        if normalize_mass:
            self.normalize_mass()
        if center_mass:
            self.center_mass()
        if center_momentum:
            self.center_momentum()
        if set_accelerations:
            self.set_accelerations()

    def do_everything(self, set_accelerations = True, normalize_mass = True, center_mass = True, center_momentum = True):
        self.do_things(set_accelerations, normalize_mass, center_mass, center_momentum)
        
    
    def add_particle(self, particle: Particle):
        self.particles.append(particle)
        for key in self.flags.keys():
            self.flags[key] = False

    #dont use _currentTime plz
<<<<<<< HEAD
    def step(self, totalTime: float):
=======
    def step_fancy(self, totalTime: float):
>>>>>>> origin/main
        smallestTime = inf
        collisionIndices = []
        collisionNumOfParticles = []
        lastParticleHadCollision = False
        for i in range(len(self.particles) - 1):
            p1 = self.particles[i]
            p2 = self.particles[i+1]
            discriminant = 2*p1.acceleration*(p2.position-p1.position) + 2*p2.acceleration*(p1.position - p2.position) + (p1.velocity-p2.velocity)**2
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
                
                if collisionTime > totalTime:       # collision doesn't happen in existing timeframe, skip
                    lastParticleHadCollision = False
                    continue
                if collisionTime > smallestTime:
                    lastParticleHadCollision = False
                    continue
                
                if collisionTime == smallestTime:        # add it to an existing collision
                    if lastParticleHadCollision:
                        collisionNumOfParticles[-1] += 1
                    else:
                        collisionIndices.append(i)
                        collisionNumOfParticles.append(2)
                        lastParticleHadCollision = True
                elif collisionTime < smallestTime:       # found new collision candidate
                    smallestTime = collisionTime
                    collisionIndices.clear()
                    collisionNumOfParticles.clear()
                    collisionIndices.append(i)
                    collisionNumOfParticles.append(2)
                    lastParticleHadCollision = True
                else:
                    continue
        
        steppingTime = min(totalTime, smallestTime)
        for particle in self.particles:
            particle.step_ghost_state(steppingTime)
        for i, collisionIndex in enumerate(collisionIndices):
            indexStart = collisionIndex
            indexEnd = collisionIndex+collisionNumOfParticles[i]
            particlesInCollision = self.particles[indexStart:indexEnd]
            newMass = sum(p.mass for p in particlesInCollision)
            newPos = particlesInCollision[0].position
            newVel = sum(p.mass*p.velocity for p in particlesInCollision)/newMass
            newParticle = Particle(newMass, newPos, newVel)
            del self.particles[indexStart:indexEnd]
            self.particles.insert(indexStart, newParticle)
            self.set_accelerations()
<<<<<<< HEAD
            self.step(totalTime-smallestTime)
=======
            self.step_fancy(totalTime-smallestTime)
>>>>>>> origin/main
            return
        return

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
              f"{"Acceleration:":13}{[round(p.acceleration, precision) for p in self.particles]}\n")
