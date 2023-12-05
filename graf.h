#pragma
#ifndef GRAF_H
#define GRAF_H

#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <fstream>
#include <unordered_set>

// creating a struct for the edges of the graph
struct edge
{
    int node1;
    int node2;
    int weight;
};

class Graph
{
private:
    int numNodes;
    std::vector<std::vector<std::pair<int, int>>> adj; // adjaency list
    std::vector<bool> visited;

    // functions for union find algorithm
    void makeSet(int node, std::vector<int> &parent, std::vector<int> &rank);
    int find(int node, std::vector<int> &parent);
    void unionSet(int node1, int node2, std::vector<int> &parent, std::vector<int> &rank);

    // function used by dfsCriticalConnections because I reset the visited vector first
    void dfsCritical(int node, int parentNode, std::vector<int> &lvl, std::vector<int> &low, std::vector<bool> &isArticulation, std::vector<std::vector<int>> &bridges);

    // resetting the visited vector
    void resetVisited() { visited = std::vector<bool>(numNodes, false); }

public:
    Graph(int numNodes = 0);

    Graph(std::vector<std::vector<std::pair<int, int>>> adj);
    Graph(int numNodes, std::vector<std::vector<int>> connections, bool isDirected = false);

    std::vector<std::vector<std::pair<int, int>>> getAdj() { return adj; }
    std::vector<bool> getVisited() { return visited; }

    // adding in the edge structure for Dijkstra algorithm
    void addEdge(int node1, int node2, int weight = 1, bool isDirected = false);

    // incresing numNodes with 1
    void addNode();

    // Depth first search returning a vector of nodes in order
    std::vector<int> dfs(int startNode);

    // DFS Critcal Connections returns a pair of vectors, the first vector contains the articulation points
    // and the second vector contains the bridges
    std::pair<std::vector<bool>, std::vector<std::vector<int>>> dfsCriticalConnections(int node);

    // Breadth first search returning a vector of nodes in order
    std::vector<int> bfs(int startNode);

    // Dijkstra's algorithm returning the shortest distance between two nodes (or two sets of nodes)
    int dijkstra(std::vector<int> startNodes, std::vector<int> endNodes);

    // Bellman Ford algorithm returning the shortest distance from one node to all the other nodes
    std::vector<int> bellmanFord(int startNode);

    // check if the graph is bipartite returning true or false
    bool isBipartite();

    // topological sorting of a graph using Kahn's algorithm returning a vector with the order of the nodes
    std::vector<int> topologicalSort();

    // minimum spanning tree using Kruskal's algorithm returning a vector of edges
    std::vector<edge> kruksalMST();

    std::vector<edge> primMST();
};

#endif // GRAF_H