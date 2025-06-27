import math
from particle import Particle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np



class ParticleSystem:
    def __init__(self, particles: list[Particle] = None, relative_velocity=0):
        if particles is None:
            self.particles = []
        else:
            self.particles = particles
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
        self.do_non_destructive()
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

    def do_non_destructive(self):
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
    
    # WARNING: THIS PERMANENTLY AFFECTS THE STATE OF THIS SYSTEM
    def step_real(self, deltaTime: float, error: float = 0.001):
        for p in self.particles:
            p.step(deltaTime)
        collisions = []
        for currentIndex in range(len(self.particles)):
            collision: set = set()
            shiftedIndex = currentIndex
            minIndex = currentIndex
            maxIndex = currentIndex
            leftBounded = False
            rightBounded = False
            for j in range(len(self.particles)):
                #top 10 code ever
                if not leftBounded:
                    leftBounded = (shiftedIndex) - 1 < 0
                    if leftBounded:
                        shiftedIndex = maxIndex
                if not rightBounded:
                    rightBounded= (shiftedIndex) + 1 >= len(self.particles)
                    if rightBounded:
                        shiftedIndex = minIndex

                if leftBounded and rightBounded:
                    continue
                elif leftBounded:
                    shiftedIndex += 1
                elif rightBounded:
                    shiftedIndex -= 1
                else:
                    shiftedIndex+= (1 - ((j & 1) << 1))*(j+1)
                    if j % 2 == 0:
                        maxIndex = shiftedIndex
                    else:
                        minIndex = shiftedIndex
                if abs(self.particles[shiftedIndex].position - self.particles[currentIndex].position) < error:
                    if len(collision) == 0:
                        collision.add(currentIndex)
                    collision.add(shiftedIndex)
                else:
                    break
            if len(collision) != 0:
                if collision not in collisions:
                    collisions.append(collision)

    if len(collisions) != 0:
        collisions = sorted(collisions, key=min)
        indexShift = 0
        for i in range(len(collisions)):
            indexStart = min(collisions[i]) + indexShift
            indexEnd = max(collisions[i]) + indexShift
            collisionParticles = system.particles[indexStart:indexEnd+1]
            mass = sum([p.mass for p in collisionParticles])
            position = sum([p.position for p in collisionParticles])/len(collisionParticles)
            velocity = sum([p.mass*p.velocity for p in collisionParticles])/mass
            for p in collisionParticles:
                system.particles.pop(indexStart)
            system.particles.insert(indexStart, Particle(mass,position,velocity))
            system.set_accelerations()
            indexShift -= indexStart + 1


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
            return [p2.mass * math.sqrt((p2.position - p1.position)/(p1.mass + p2.mass)), (-1) * (p1.mass) * math.sqrt((p2.position - p1.position)/(p1.mass + p2.mass))]
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
