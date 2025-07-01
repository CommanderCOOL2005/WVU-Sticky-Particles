from random import *

from particle import *
from particle_system import *
from plotter import *

innputte = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(4)]
system = ParticleSystem(innputte)
system.assign_total_question_mark_solution()
plot(system, 3, 200)
