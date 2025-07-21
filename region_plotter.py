import colorsys
import csv
import sys

import pygame
from pygame.locals import *

from plotter import *

print("File to read: ", end="")
filename = input()
solution = None
TARGET_INDEX = None
IMAGE_WIDTH = None
IMAGE_HEIGHT = None
POS_MIN_SHIFT = None
POS_MAX_SHIFT = None
VEL_MIN_SHIFT = None
VEL_MAX_SHIFT = None
data = []
with open(file=f"region_data/{filename}.csv", mode='r') as csvfile:
    masses = None
    positions = None
    velocities = None
    reader = csv.reader(csvfile)
    for row in reader:
        if reader.line_num == 1:
            masses = [float(elem) for elem in row]
        elif reader.line_num == 2:
            positions = [float(elem) for elem in row]
        elif reader.line_num == 3:
            velocities = [float(elem) for elem in row]
        elif reader.line_num == 4:
            TARGET_INDEX = int(row[0])
            IMAGE_WIDTH = int(row[1])
            IMAGE_HEIGHT = int(row[2])
            POS_MIN_SHIFT = float(row[3])
            POS_MAX_SHIFT = float(row[4])
            VEL_MIN_SHIFT = float(row[5])
            VEL_MAX_SHIFT = float(row[6])
        else:
            data.append([int(row[0]), int(row[1]), True if row[2] == "True" else False, float(row[3])])
    initial_input = []
    for i in range(len(masses)):
        initial_input.append(Particle(masses[i],positions[i],velocities[i]))
    solution = ParticleSystem(initial_input)

pygame.init()
 
fps = 60
fps_clock = pygame.time.Clock()
 
width, height = IMAGE_WIDTH, IMAGE_HEIGHT
screen = pygame.display.set_mode((width, height))
mouse_pos = None
mouse_last_click = None
mouse_is_clicking = False

region_surface = pygame.Surface((IMAGE_WIDTH, IMAGE_HEIGHT))
#Precompute the surface
min_ham = min(p[3] for p in data)
max_ham = max(p[3] for p in data)
for p in data:
    color = (0,0,0)
    heat = (p[3] - min_ham)/(max_ham - min_ham)
    if p[2]:
        color = [round(255*v) for v in colorsys.hsv_to_rgb(0.35 - 0.2*heat,0.8*(1-heat),0.9)]
    else:
        color = [120*heat, 120*heat, 120*heat]
    region_surface.set_at((p[0],p[1]),color)

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
        pygame.draw.circle(screen, center=mouse_pos, radius=4, color=(255,255,255))
    
    if mouse_is_clicking and mouse_last_click:
        world_space_vector = (
            map_range(mouse_last_click[0], 0, IMAGE_WIDTH, POS_MIN_SHIFT, POS_MAX_SHIFT),
            map_range(mouse_last_click[1], IMAGE_HEIGHT, 0, VEL_MIN_SHIFT, VEL_MAX_SHIFT)
        )
        sol_to_plot = copy.deepcopy(solution)
        sol_to_plot.shift_particle(TARGET_INDEX, position=world_space_vector[0], velocity=world_space_vector[1])
        plot_solution(sol_to_plot, 3, 200)

    pygame.display.flip()  
    fps_clock.tick(fps)