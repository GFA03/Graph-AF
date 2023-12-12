from build.graf import *
import pygame
import sys

graph = Graph(5)

graph.addEdge(0, 1, 4, False)
graph.addEdge(0, 2, 3, False)
graph.addEdge(0, 3, 2, False)
graph.addEdge(0, 4, 1, False)
graph.addEdge(1, 2, 1, False)
graph.addEdge(1, 3, 2, False)
graph.addEdge(1, 4, 3, False)
graph.addEdge(2, 3, 4, False)
graph.addEdge(2, 4, 5, False)
graph.addEdge(3, 4, 6, False)

pygame.init()

screen = pygame.display.set_mode((1280, 720))
pygame.display.set_caption("Graphs")
clock = pygame.time.Clock()
running = True
dt = 0

player_pos = pygame.Vector2(screen.get_width() / 2, screen.get_height() / 2)

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    
    screen.fill('lightgray')

    pygame.draw.circle(screen, "red", player_pos, 40)

    keys = pygame.key.get_pressed()
    if keys[pygame.K_w]:
        player_pos.y -= 300 * dt
    if keys[pygame.K_s]:
        player_pos.y += 300 * dt
    if keys[pygame.K_a]:
        player_pos.x -= 300 * dt
    if keys[pygame.K_d]:
        player_pos.x += 300 * dt

    pygame.display.flip()

    dt = clock.tick(60) / 1000

pygame.display.quit()
pygame.quit()