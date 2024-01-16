# Graph Algorithms Project

This project is an implementation of various graph algorithms in C++. The goal is to provide an efficient and versatile set of tools for working with graphs. Currently, the project includes implementations of the following graph algorithms:

## Implemented Graph Algorithms

### Depth-First Search (DFS)

The Depth-First Search algorithm is utilized for traversing or searching graph data structures. It explores as far as possible along each branch before backtracking.

### Breadth-First Search (BFS)

Breadth-First Search is employed for traversing or searching graph data structures. It systematically explores the neighbor nodes at the present depth before moving on to nodes at the next depth level.

### DFS for Finding Critical Nodes (Articulation Points)

This algorithm identifies critical nodes, also known as articulation points, in a graph. Articulation points are those whose removal would increase the number of connected components in the graph.

### Dijkstra's Algorithm

Dijkstra's Algorithm is a popular algorithm for finding the shortest path between nodes in a graph. It is particularly useful in scenarios where weighted edges are present.

### Bellman-Ford Algorithm

The Bellman-Ford Algorithm is employed for finding the shortest paths from a single source vertex to all other vertices, even in the presence of negative edge weights.

### Check if a Graph is Bipartite

This algorithm determines whether a given graph is bipartite, meaning it can be colored using two colors in such a way that no two adjacent nodes have the same color.

### Finding an Eulerian Path

The Eulerian Path algorithm identifies a path in a graph that visits every edge exactly once. If such a path exists, the graph is considered Eulerian.

## License

This project is licensed under the [MIT License](LICENSE).
