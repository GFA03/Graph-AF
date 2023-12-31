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
#include <unordered_map>

#define INT_MAX 1e9

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
    int out[1000];

    // functions for union find algorithm
    void makeSet(int node, std::vector<int> &parent, std::vector<int> &rank);
    int find(int node, std::vector<int> &parent);
    void unionSet(int node1, int node2, std::vector<int> &parent, std::vector<int> &rank);

    // function used by dfsCriticalConnections because I reset the visited vector first
    void dfsCritical(int node, int parentNode, std::vector<int> &lvl, std::vector<int> &low, std::vector<bool> &isArticulation, std::vector<std::vector<int>> &bridges);

    // resetting the visited vector
    void resetVisited() { visited = std::vector<bool>(numNodes, false); }

    // function for counting the number of in and out degrees and the number of edges
    void countInOutDegrees(int in[], int out[], int &edges);

    // function for finding the starting node of an Eulerian Path
    int findStartNode();

public:
    Graph(int numNodes = 0);

    Graph(std::vector<std::vector<int>> &adj);
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

    // DFS made for finding the Eulerian path in a graph
    void dfsEuler(int startNode, int outTemp[], std::deque<int> &path);

    // Breadth first search returning a vector of nodes in order
    std::vector<int> bfs(int startNode);

    // Breadth first search that allows repetition of nodes and avoids indefinite loops by checking if the state was already visited
    int bfsBitManipulation();

    // Dijkstra's algorithm returning the shortest distance between two nodes (or two sets of nodes)
    int dijkstra(std::vector<int> startNodes, std::vector<int> endNodes);

    // Bellman Ford algorithm returning the shortest distance from one node to all the other nodes
    std::vector<int> bellmanFord(int startNode);

    // check if the graph is bipartite returning true or false
    bool isBipartite();

    // topological sorting of a graph using Kahn's algorithm returning a vector with the order of the nodes
    std::vector<int> topologicalSort();

    // minimum spanning tree using Kruskal's algorithm returning a vector of edges
    std::vector<edge> kruskalMST();

    std::vector<edge> primMST();

    bool isEulerian();

    std::deque<int> findEulerianPath();
};

#endif // GRAF_H