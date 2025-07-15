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
                _relative_velocity=0
                ):
        """Creates a system of particles that represent when the initial mass measure (rho_0) in the Repulsive Pressureless Euler Poisson system is a weighted sum of dirac masses.
        
        The inputted particles are automatically adjusted so that:
            1) Their total mass is 1. 
            2) Their center of mass is 0.
            3) Their total momentum is 0,
            4) Their accelerations are set according to RPEP.
        If for some reason you do not want these adjustments, you can override them by passing in override_normalization = True (which will override 1,2,3) or override_acceleration = True (which will override 4). You can configure these adjustments with the methods configure(), normalize_system() and set_accelerations()

        Args:
            particles (list[Particle]): A list of particles 
            override_normalization (bool, optional): Whether to override setting total mass to 1, center of mass to 0, and total momentum to 0. Defaults to False.
            override_acceleration (bool, optional): Whether to override setting their accelerations according to RPEP. Defaults to False.
        """
        if particles is None:
            self.particles = []
        else:
            self.particles = sorted(particles, key=lambda particle: particle.position)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
        self.make_rgb_colors()
        self.rel_vel = _relative_velocity  # relative velocity of the system

    def get_total_mass(self):
        return sum(p.mass for p in self.particles)
    
    def get_center_of_mass(self):
        return sum(p.mass * p.position for p in self.particles) / self.get_total_mass()
        
    def get_total_momentum(self):
        return sum(p.mass * p.velocity for p in self.particles)

    def set_accelerations(self):
        """Sets the acceleration of every particle so that they follow the RPEP system
        """
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
        """Adjusts the particles in the system so that their sum is 1.
        """
        total_mass = self.get_total_mass()
        for p in self.particles:
            p.mass /= total_mass

    def center_mass(self):
        """Adjusts the particles in the system so that their center of mass is 0.
        """
        com = self.get_center_of_mass()
        for p in self.particles:
            p.position -= com

    def center_momentum(self):
        """Adjusts the particles in the system so that their center of mass is 0.
        """
        total_momentum = sum(p.mass * p.velocity for p in self.particles)
        total_mass = sum(p.mass for p in self.particles)
        mean_velocity = total_momentum / total_mass

        for p in self.particles:
            p.velocity -= mean_velocity
    
    def configure(self, normalize_mass = False, center_mass = False, center_momentum = False):
        """Configures the system in various ways

        Args:
            normalize_mass (bool, optional): Calls :py:func:`normalize_mass()` if true. Defaults to False.
            center_mass (bool, optional): Calls :py:func:`center_mass()` if true. Defaults to False.
            center_momentum (bool, optional): Calls :py:func:`center_momentum()` if true. Defaults to False.
        """
        self.particles = sorted(self.particles, key=lambda particle: particle.position)
        if normalize_mass:
            self.normalize_mass()
        if center_mass:
            self.center_mass()
        if center_momentum:
            self.center_momentum()

    def normalize_system(self, normalize_mass = True, center_mass = True, center_momentum = True):
        """Normalizes the system so that the system's total mass is one, center of mass is zero, and total momentum is zero.

        Args:
            normalize_mass (bool, optional): Defaults to True.
            center_mass (bool, optional): Defaults to True.
            center_momentum (bool, optional): Defaults to True.
        """
        self.configure(normalize_mass, center_mass, center_momentum)
    
    def add_particle(self, particle: Particle, override_normalization = False, override_acceleration = False):
        """Adds a new particle to the system and renormalizes the system and assigns new accelerations to each particle

        Args:
            particle (Particle): _description_
            override_normalization (bool, optional): Defaults to False.
            override_acceleration (bool, optional): Defaults to False.
        """
        self.particles.append(particle)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
    
    def add_particles(self, particles: list[Particle], override_normalization = False, override_acceleration = False):
        """Adds multiple new particles to the system and renormalizes the system.

        Args:
            particles (list[Particle]): _description_
            override_normalization (bool, optional):  Defaults to False.
            override_acceleration (bool, optional): Defaults to False.
        """
        self.particles.extend(particles)
        if not override_normalization:
            self.normalize_system()
        if not override_acceleration:
            self.set_accelerations()
    
    def get_next_collision(self) -> Collision:
        """Gets the next collision that the system will experience. The information about the collision is stored in a ``Collision`` object. 

        Returns:
            Collision: The next collision the system will experience. If the system experiences no collision ever again, this method returns a ``Collision`` object with infinite time.
        """
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
                
                error = 1e-15
                if abs(collisionTime-earliest.time) < error:        # add it to an existing collision
                    if lastParticleHadCollision:
                        earliest.sizes[-1] += 1
                    else:
                        earliest.indices.append(i)
                        earliest.sizes.append(2)
                        lastParticleHadCollision = True
                elif collisionTime + leeway >= earliest.time:
                    lastParticleHadCollision = False
                    continue
                elif 0 < collisionTime - leeway <= earliest.time:       # found new collision candidate, exclude zero
                    earliest.time = collisionTime
                    earliest.indices.clear()
                    earliest.sizes.clear()
                    earliest.indices.append(i)
                    earliest.sizes.append(2)
                    lastParticleHadCollision = True
                else:
                    continue
        return earliest

    def advance(self, total_time: float, _collision : Collision = None):
        """Simulates the system for ``total_time`` time. This *permanently* affects the state of this system. 

        Args:
            total_time (float): The amount of time to simulate the system for
        """
        if total_time <= 0:
            return
        if _collision == None:
            collision= self.get_next_collision()
        else:
            collision = _collision
        
        steppingTime = min(total_time, collision.time)
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

            self.advance(total_time-collision.time)
            return
        return

    def advance_to_next_collision(self) -> float:
        """Simulates the system up until the next collision. This *permanently* affects the state of this system.

        Returns:
            float: The time of collision
        """
        collision = self.get_next_collision()
        time = 0
        if collision.time != inf:
            time = collision.time
        self.advance(time, collision)
        return time
    
    def is_equilibrium_solution(self) -> bool:
        """Checks if this system converges to equilibrium in finite time

        Returns:
            bool: Whether or not it converges to equilibrium in finite time
        """
        system_copy = copy.deepcopy(self)
        while True:
            time = system_copy.advance_to_next_collision()
            if time == 0:
                if len(system_copy.particles) > 1:
                    return False
                else:
                    return True
        

    
    def next_collision(self):
        # Find the time of the next perfect collision in the system.
        # Each collision gives us the following information:
        # v1 - v2 = sqrt((y2 - y1) * (m1 + m2))
        # After the collision, the particles will combined into one larger particle.
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
        """Sets the velocity of every particle in the system so that the solution becomes perfect
        """
        velocities = self.get_perfect_solution()
        for i in range(len(self.particles)):
            self.particles[i].velocity = velocities[i]
        self.center_momentum()
    
    
    def assign_random_signed_velocities(self, a: float = 0, b:float = 1):
        """Sets the magnitude of the velocity of each particle between ``a`` and ``b``. The sign of the velocity will be positive if the position of the particle is less than zero, and be negative if the position of the particle is greater than zero 

        Args:
            a (float, optional): _description_. Defaults to 0.
            b (float, optional): _description_. Defaults to 1.
        """
        for i in range(len(self.particles)):
            self.particles[i].velocity = (1 if self.particles[i].position < 0 else -1)*uniform(a, b)
        self.center_momentum()

    def make_rgb_colors(self):
        """Makes it so the particles have a rainbow color pallete when graphed.
        """
        for i, particle in enumerate(self.particles):
            particle.color = matplotlib.colors.hsv_to_rgb((i/len(self.particles),1,0.75))

    def print_info(self):
        """Prints some characteristic and debug information about the system
        """
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
        """Prints information about each particle
        """
        precision = 5
        print(f"{"Mass:":13}{[round(p.mass, precision) for p in self.particles]}\n" +
              f"{"Position:":13}{[round(p.position, precision) for p in self.particles]}\n"+
              f"{"Velocity:":13}{[round(p.velocity, precision) for p in self.particles]}\n"+
              f"{"Acceleration:":13}{[round(p.acceleration, precision) for p in self.particles]}")
    
    def dump(self):
        for i in range(len(self.particles)):
            print(self.particles[i], end=',')
        print(self.particles[-1])