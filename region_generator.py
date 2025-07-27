import copy
import csv
import os
 
from plotter import *

while True:
    initial_input = [Particle(uniform(0.01,1), uniform(-1,1),2) for i in range(20)]
    solution = ParticleSystem(initial_input)
    solution.assign_perfect_solution()
    for p in solution.particles:
        p.velocity *= uniform(1,1.2)
    solution.adjust_solution()
    plot_solution(system=solution, plot_until_time_bound=True)
    confirmation = input()
    if confirmation == "y":
        break


TARGET_INDEX = 0

IMAGE_WIDTH = 1000
IMAGE_HEIGHT = 1000
POS_MIN_SHIFT = -2
POS_MAX_SHIFT = 2
POS_DELTA = (POS_MAX_SHIFT - POS_MIN_SHIFT)/IMAGE_WIDTH
VEL_MIN_SHIFT = -2
VEL_MAX_SHIFT = 2
VEL_DELTA = (VEL_MAX_SHIFT - VEL_MIN_SHIFT)/IMAGE_HEIGHT

data : list[dict] = []

for i in range(IMAGE_HEIGHT):
    print(f"{i+1}/{IMAGE_HEIGHT}")
    for j in range(IMAGE_WIDTH):
        shifted_sol = copy.deepcopy(solution)
        shifted_sol.shift_particle(
           particle_index=0, 
           position=POS_MIN_SHIFT + POS_DELTA*j,
           velocity=VEL_MIN_SHIFT + VEL_DELTA*i
           )
        data.append({
                "x": j,
                "y": IMAGE_HEIGHT - i,
                "eq": shifted_sol.is_equilibrium_solution(),
                "ham": shifted_sol.get_hamiltonian() - shifted_sol.get_total_momentum()**2
            })

if not os.path.exists("region_data"):
    os.makedirs("region_data")

filename = f"region_data/{hash(solution)}.csv"

with open(filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([p.mass for p in solution.particles])
    writer.writerow([p.position for p in solution.particles])
    writer.writerow([p.velocity for p in solution.particles])
    writer.writerow([TARGET_INDEX, IMAGE_WIDTH, IMAGE_HEIGHT, POS_MIN_SHIFT, POS_MAX_SHIFT, VEL_MIN_SHIFT, VEL_MAX_SHIFT])
    for elem in data:
        writer.writerow(elem.values())