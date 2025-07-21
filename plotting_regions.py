import copy
import sys

import matplotlib.pyplot as plt
import pygame
from pygame.locals import *
 
from particle import *
from particle_system import *
from plotter import *

def map_range(x, input_start, input_end, output_start, output_end):
    return  (x - input_start) / (input_end - input_start) * (output_end - output_start) + output_start


initial_input = [Particle(uniform(0.01,1), uniform(-1,1)) for i in range(20)]
# initial_input = [Particle(0.5,-1), Particle(0.5,1)]
solution = ParticleSystem(initial_input)
solution.assign_perfect_solution()
for p in solution.particles:
    p.velocity *= 1.00000001
solution.adjust_solution()
solution.print_particles()

target_index = 0

IMAGE_WIDTH = 800
IMAGE_HEIGHT = 800
pos_min_shift = -1
pos_max_shift = 2
pos_range = pos_max_shift - pos_min_shift
pos_delta = pos_range/IMAGE_WIDTH
vel_min_shift = -3
vel_max_shift = 3
vel_range = vel_max_shift - vel_min_shift
vel_delta = vel_range/IMAGE_HEIGHT

big_bertha = []

for i in range(IMAGE_HEIGHT):
    print(f"{i+1}/{IMAGE_HEIGHT}")
    for j in range(IMAGE_WIDTH):
        shifted_sol = copy.deepcopy(solution)
        shifted_sol.shift_particle(
           particle_index=0, 
           position=pos_min_shift + pos_delta*j,
           velocity=vel_min_shift + vel_delta*i
           )
        big_bertha.append({
                "x": j,
                "y": IMAGE_HEIGHT - i,
                "eq": shifted_sol.is_equilibrium_solution(),
                "ham": shifted_sol.get_hamiltonian() - shifted_sol.get_total_momentum()**2
            })

# print(big_bertha)

pygame.init()
 
fps = 60
fpsClock = pygame.time.Clock()
 
width, height = IMAGE_WIDTH, IMAGE_HEIGHT
screen = pygame.display.set_mode((width, height))
mouse_pos = None
mouse_last_click = None
mouse_is_clicking = False

region_surface = pygame.Surface((IMAGE_WIDTH, IMAGE_HEIGHT))
#Precompute the surface
min_ham = min(p["ham"] for p in big_bertha)
max_ham = max(p["ham"] for p in big_bertha)
for p in big_bertha:
    color = (0,0,0)
    if p["eq"]:
        heat = 255*(p["ham"] - min_ham)/(max_ham - min_ham)
        color = (heat, 100, 255-heat)
    region_surface.set_at((p["x"],p["y"]),color)

# Game loop.
while True:
    screen.fill((0, 0, 0))
  
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()
        if event.type == MOUSEMOTION:
            mouse_pos = pygame.mouse.get_pos()
        if event.type == MOUSEBUTTONUP:
            mouse_last_click = pygame.mouse.get_pos()
            mouse_is_clicking = True
        else:
            mouse_is_clicking = False
    
    screen.blit(region_surface, (0,0))
    if mouse_pos:
        pygame.draw.circle(screen, center=mouse_pos, radius=2, color=(255,255,255))
    
    if mouse_is_clicking and mouse_last_click:
        world_space_vector = (
            map_range(mouse_last_click[0], 0, IMAGE_WIDTH, pos_min_shift, pos_max_shift),
            map_range(mouse_last_click[1], IMAGE_HEIGHT, 0, vel_min_shift, vel_max_shift)
        )
        sol_to_plot = copy.deepcopy(solution)
        sol_to_plot.shift_particle(target_index, position=world_space_vector[0], velocity=world_space_vector[1])
        plot_solution(sol_to_plot, 3, 200)

    pygame.display.flip()  
    fpsClock.tick(fps)