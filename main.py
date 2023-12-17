from build.graf import *
import pygame
import sys

class GraphVisualization:
    def __init__(self, graph):
        self.graph = graph
        self.node_radius = 20
        self.node_color = (255, 0, 0)  # RGB color for nodes
        self.node_text_color = (255, 255, 255)  # RGB color for node text
        self.edge_color = (0, 0, 0)   # RGB color for edges
        self.message_color = (255, 0, 0)  # RGB color for alert message

        self.add_node_button_rect = pygame.Rect(10, 10, 100, 30)
        self.add_node_button_color = (0, 255, 0)  # RGB color for the button
        self.add_node_button_text_color = (255, 255, 255)  # RGB color for the button text

        self.adding_node = False
        self.collision_alert = False

        self.dragging = False
        self.selected_node_index = None
        self.offset_x, self.offset_y = 0, 0

    def does_node_intersect(self, pos, exclude_index = None):
        for node, coordinates in self.graph.getCoordinates().items():
            if node == exclude_index:
                continue
            distance = ((pos[0] - coordinates[0]) ** 2 + (pos[1] - coordinates[1]) ** 2) ** 0.5
            if distance < 2 * self.node_radius:  # Assuming the nodes are circles with radius node_radius
                return True
        return False

    def draw_graph(self):
        pygame.init()

        # Set up the display
        width, height = 800, 600
        screen = pygame.display.set_mode((width, height))
        pygame.display.set_caption("Graph Visualization")

        # Main loop
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False

                elif event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == 1:  # Left mouse button
                        # Check if the add node button is clicked
                        if self.add_node_button_rect.collidepoint(event.pos):
                            self.adding_node = True
                        elif self.adding_node:
                            if not self.does_node_intersect(event.pos):
                                self.graph.addNode([event.pos[0], event.pos[1]])
                                self.adding_node = False
                            else:
                                # Display collision alert
                                self.collision_alert = True

                        else:
                            for i, node_coordinates in enumerate(self.graph.getCoordinates().values()):
                                if self.is_point_inside_circle(event.pos, node_coordinates, self.node_radius):
                                    self.dragging = True
                                    self.selected_node_index = len(self.graph.getAdj()) - i - 1
                                    self.offset_x = node_coordinates[0] - event.pos[0]
                                    self.offset_y = node_coordinates[1] - event.pos[1]

                elif event.type == pygame.MOUSEBUTTONUP:
                    if event.button == 1:
                        self.dragging = False
                        self.selected_node_index = None

                elif event.type == pygame.MOUSEMOTION:
                    if self.dragging and self.selected_node_index is not None:
                        new_x = event.pos[0] + self.offset_x
                        new_y = event.pos[1] + self.offset_y

                        if not self.does_node_intersect((new_x, new_y), exclude_index = self.selected_node_index):
                            self.graph.setCoordinates(self.selected_node_index, (new_x, new_y))

            screen.fill((255, 255, 255))  # Fill the background with white

            # Draw edges
            # (Draw edges logic)

            # Draw nodes
            for node, coordinates in self.graph.getCoordinates().items():
                pygame.draw.circle(screen, self.node_color, (int(coordinates[0]), int(coordinates[1])), self.node_radius)
                font = pygame.font.Font(None, 24)
                text = font.render(str(node), True, self.node_text_color)
                text_rect = text.get_rect(center=(int(coordinates[0]), int(coordinates[1])))
                screen.blit(text, text_rect)

            # Draw the add node button
            pygame.draw.rect(screen, self.add_node_button_color, self.add_node_button_rect)
            font = pygame.font.Font(None, 24)
            text = font.render("Add Node", True, self.add_node_button_text_color)
            screen.blit(text, (self.add_node_button_rect.x + 10, self.add_node_button_rect.y + 5))

            # Display collision alert
            if self.collision_alert:
                alert_font = pygame.font.Font(None, 24)
                alert_text = alert_font.render("Collision! Choose another position.", True, self.message_color)
                screen.blit(alert_text, (width // 2 - alert_text.get_width() // 2, height - 50))

            pygame.display.flip()

        pygame.quit()
        sys.exit()

    def is_point_inside_circle(self, point, circle, radius):
        return (point[0] - circle[0]) ** 2 + (point[1] - circle[1]) ** 2 <= radius ** 2


# Example usage
graph = Graph(0)
visualization = GraphVisualization(graph)
visualization.draw_graph()