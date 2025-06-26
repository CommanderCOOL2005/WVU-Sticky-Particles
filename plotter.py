import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from particle_system import *

def get_x_bounds(system: ParticleSystem): 
        positions = [p.position for p in system.particles]
        bound = max(abs(min(positions)),abs(max(positions)))
        bound *= 1.5
        return (-bound, bound)

def plot_ghost_state(system: ParticleSystem):
    height = 15
    fig, ax = plt.subplots()
    ax.set(
        xlim=(get_x_bounds(system)),
        ylim=(0,height), 
        xlabel=r'Position'
    )

    for i in range(len(system.particles)):
        y = np.linspace(0,height,100)
        x = [system.particles[i].evaluate_ghost_state(time) for time in y]
        ax.plot(x,y, color=matplotlib.colors.hsv_to_rgb((i/len(system.particles),1,1)), lw=0.5)
    plt.show()

def plot_real_state(system: ParticleSystem, end_time: float, steps: int):
    deltaTime = end_time/steps
    fig, ax = plt.subplots()
    ax.set(
        xlim=(get_x_bounds(system)),
        ylim=(0,end_time+1)
    )
    particle_histories = {} 
    for p in system.particles:
        particle_histories[p.id] = [(0, p.position)]
    for i in range(steps):
        step_real(system, deltaTime)
        current_time = i * deltaTime
        for p in system.particles:
            if p.id not in particle_histories:
                particle_histories[p.id] = []
            particle_histories[p.id].append((current_time,p.position))
    for history in particle_histories.values():
        times, positions = zip(*history)
        ax.plot(positions, times, lw=0.5)
    plt.show()

# WARNING: THIS PERMANENTLY AFFECTS THE STATE OF THIS SYSTEM
def step_real(system: ParticleSystem, deltaTime: float, error: float = 0.001):
    for p in system.particles:
        p.velocity += p.acceleration * deltaTime
        p.position += p.velocity * deltaTime
    collisions = []
    for currentIndex in range(len(system.particles)):
        collision: set = set()
        shiftedIndex = currentIndex
        minIndex = currentIndex
        maxIndex = currentIndex
        leftBounded = False
        rightBounded = False
        for j in range(len(system.particles)):
            #top 10 code ever
            if not leftBounded:
                leftBounded = (shiftedIndex) - 1 < 0
                if leftBounded:
                    shiftedIndex = maxIndex
            if not rightBounded:
                rightBounded= (shiftedIndex) + 1 >= len(system.particles)
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
            if abs(system.particles[shiftedIndex].position - system.particles[currentIndex].position) < error:
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