#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <fstream>

// creating a struct for the edges of the graph
struct edge {
    int node1;
    int node2;
    int weight;
};

// overloading the < operator for the edges
bool operator < (const edge& a,const edge& b) {
    return a.weight < b.weight;
}

class Graph {
private:
    int numNodes;
    std::vector<std::vector<int>> adj;
    std::vector<edge> edges;
    std::vector<bool> visited;
    std::vector<int> lvl;
    std::vector<int> low;
    std::vector<bool> isArticulation;
    int dfsCriticalConditions(int node, int parentNode);
    void resetVisited();
    void resetLvl();
    void resetLow();
    void resetIsArticulation();
    void resetAll();
public:
    Graph(int numNodes = 0);
    Graph(int numNodes, std::vector<std::vector<int>> adj);
    Graph(int numNodes, std::vector<std::vector<int>> adj, std::vector<edge> edges);

    void addEdge(int node1, int node2, int weight = 1, bool isDirected = false);
    void addNode();
    std::vector<int> dfs(int startNode);
    std::vector<int> bfs(int startNode);
    int dijkstra(int startNode, int endNode, std::vector<edge> edges);
    int dijkstra(int startNode, int endNode);
    std::vector<int> bellmanFord(int startNode);
    bool isBipartite(int node, int parentNode);
    std::vector<int> topologicalSort();
    std::vector<int> minimumSpanningTree();
};
// initialising graph using just the number of nodes
Graph::Graph(int numNodes) {
    this->numNodes = numNodes;
    adj.resize(numNodes);
    visited.resize(numNodes, false);
    lvl.resize(numNodes, -1);
    low.resize(numNodes, -1);
    isArticulation.resize(numNodes, false);
}

// initialising graph using the number of nodes and the adjacency list
Graph::Graph(int numNodes, std::vector<std::vector<int>> adj) {
    this->numNodes = numNodes;
    this->adj = adj;
    // for(int i = 0; i < adj.size(); ++i) {
    //     for(int j = 0; j < adj[i].size(); ++j) {
    //         edges.push_back({i, adj[i][j], 1});
    //     }
    // }
    visited.resize(numNodes, false);
    lvl.resize(numNodes, -1);
    low.resize(numNodes, -1);
    isArticulation.resize(numNodes, false);
}   

// initialising graph using the number of nodes, the adjacency list and the vector of edges
Graph::Graph(int numNodes, std::vector<std::vector<int>> adj, std::vector<edge> edges){
    this->numNodes = numNodes;
    this->adj = adj;
    this->edges = edges;
    visited.resize(numNodes, false);
    lvl.resize(numNodes, -1);
    low.resize(numNodes, -1);
    isArticulation.resize(numNodes, false);
}

void Graph::addNode(){
    numNodes++;
}

// adding edge to the graph
void Graph::addEdge(int node1, int node2, int weight, bool isDirected) {
    if(!isDirected) {
        adj[node1].push_back(node2);
        adj[node2].push_back(node1);
        edges.push_back({node1, node2, weight});
        return;
    }
    adj[node1].push_back(node2);
    edges.push_back({node1, node2, weight});
}

// reseting the visited vector
void Graph::resetVisited() {
    for (int i = 0; i < numNodes; i++) {
        visited[i] = false;
    }
}

// reseting the lvl vector
void Graph::resetLvl() {
    for (int i = 0; i < numNodes; i++) {
        lvl[i] = -1;
    }
}

// reseting the low vector
void Graph::resetLow() {
    for (int i = 0; i < numNodes; i++) {
        low[i] = -1;
    }
}

// reseting the isArticulation vector
void Graph::resetIsArticulation() {
    for (int i = 0; i < numNodes; i++) {
        isArticulation[i] = false;
    }
}

// reseting all vectors
void Graph::resetAll() {
    resetVisited();
    resetLvl();
    resetLow();
    resetIsArticulation();
}

// Depth first search returning a vector with the order of the nodes visited
std::vector<int> Graph::dfs(int startNode) {
    resetAll();
    std::stack<int> s;
    std::vector<int> returningDFS;
    visited[startNode] = true;
    s.push(startNode);
    returningDFS.push_back(startNode);
    while(!s.empty()) {
        int currentNode = s.top();
        s.pop();
        returningDFS.push_back(currentNode);
        for (int i = 0; i < adj[currentNode].size(); i++) {
            if (!visited[adj[currentNode][i]]) {
                visited[adj[currentNode][i]] = true;
                s.push(adj[currentNode][i]);
            }
        }
    }
    return returningDFS;
}

// Breadth first search returning a vector with the order of the nodes visited
std::vector<int> Graph::bfs(int startNode) {
    resetAll();
    std::vector<int> returningBFS;
    std::queue<int> q;
    visited[startNode] = true;
    q.push(startNode);
    while(!q.empty()) {
        int currentNode = q.front();
        q.pop();
        returningBFS.push_back(currentNode);
        for (int i = 0; i < adj[currentNode].size(); i++) {
            if (!visited[adj[currentNode][i]]) {
                visited[adj[currentNode][i]] = true;
                q.push(adj[currentNode][i]);
            }
        }
    }
    return returningBFS;
}

// TO DO: Dijkstra algorithm should return VECTOR not INT and you can get the SIZE by VECTOR.SIZE() instead
// Dijkstra's algorithm returning the shortest distance between two nodes using a vector of edges
int Graph::dijkstra(int startNode, int endNode) {
    std::vector<int> dist(numNodes, 1000000000);
    std::priority_queue<std::pair<int, int>> pq;
    if(edges.empty())
        for(int i = 0; i < adj.size(); ++i) {
            for(int j = 0; j < adj[i].size(); ++j) {
                edges.push_back({i, adj[i][j], 1});
            }
    }
    dist[startNode] = 0;
    pq.push(std::make_pair(0, startNode));
    while(!pq.empty()) {
        int currentNode = pq.top().second;
        pq.pop();
        for (int i = 0; i < edges.size(); i++) {
            if (edges[i].node1 == currentNode) {
                int nextNode = edges[i].node2;
                int weight = edges[i].weight;
                if (dist[nextNode] > dist[currentNode] + weight) {
                    dist[nextNode] = dist[currentNode] + weight;
                    pq.push(std::make_pair(-dist[nextNode], nextNode));
                }
            }
        }
    }
    return dist[endNode];
}

std::vector<int> Graph::bellmanFord(int startNode) {
    std::vector<int> dist(numNodes, 1000000000);
    dist[startNode] = 0;
    for (int i = 0; i < numNodes - 1; i++) {
        for (int j = 0; j < edges.size(); j++) {
            int node1 = edges[j].node1;
            int node2 = edges[j].node2;
            int weight = edges[j].weight;
            if (dist[node2] > dist[node1] + weight) {
                dist[node2] = dist[node1] + weight;
            }
        }
    }
    return dist;
}

// check if the graph is bipartite, before using it you should call ResetAll() method because the visited vector might be changed
bool Graph::isBipartite(int node, int parentNode) {
    std::queue<int> q;
    q.push(node);
    std::vector<int> colors(numNodes, -1);
    colors[node] = 0;
    while(!q.empty()) {
        int currentNode = q.front();
        q.pop();
        for(int i = 0; i < adj[currentNode].size(); ++i) {
            int nextNode = adj[currentNode][i];
            if(colors[nextNode] == -1) {
                colors[nextNode] = 1 - colors[currentNode];
                q.push(nextNode);
            } else if(colors[nextNode] == colors[currentNode]) {
                return false;
            }
        }
    }
    return true;
}

// topological sorting of a graph using Kahn's algorithm returning a vector with the order of the nodes
std::vector<int> Graph::topologicalSort() {
    std::vector<int> returningTopologicalSort;
    std::vector<int> inDegree(numNodes, 0);
    for(int i = 0; i < numNodes; ++i) {
        for(int j = 0; j < adj[i].size(); ++j) {
            inDegree[adj[i][j]]++;
        }
    }
    std::queue<int> q;
    for(int i = 0; i < numNodes; ++i) {
        if(inDegree[i] == 0) {
            q.push(i);
        }
    }
    while(!q.empty()) {
        int currentNode = q.front();
        q.pop();
        returningTopologicalSort.push_back(currentNode);
        for(int i = 0; i < adj[currentNode].size(); ++i) {
            int nextNode = adj[currentNode][i];
            inDegree[nextNode]--;
            if(inDegree[nextNode] == 0) {
                q.push(nextNode);
            }
        }
    }
    return returningTopologicalSort;
}

// minimum spanning tree
std::vector<int> Graph::minimumSpanningTree() {
    std::vector<int> returningMST;
    std::vector<edge> sortedEdges = edges;
    // check if sortedEdges does a good job in sorting the edges
    std::sort(sortedEdges.begin(), sortedEdges.end());
    std::vector<int> parent(numNodes);
    for(int i = 0; i < numNodes; ++i) {
        parent[i] = i;
    }
    for(int i = 0; i < sortedEdges.size(); ++i) {
        int node1 = sortedEdges[i].node1;
        int node2 = sortedEdges[i].node2;
        int weight = sortedEdges[i].weight;
        if(parent[node1] != parent[node2]) {
            returningMST.push_back(weight);
            int oldParent = parent[node1];
            int newParent = parent[node2];
            for(int j = 0; j < numNodes; ++j) {
                if(parent[j] == oldParent) {
                    parent[j] = newParent;
                }
            }
        }
    }
    return returningMST;
}

// articulation points and bridges

// deseneaza pe o foaie Graphuri ca sa intelegi articulation points and bridges
    
int main () {
    std::ifstream fin("Graph.in");
    int n, m;
    fin >> n >> m;
    std::vector<std::vector<int>> adj(n+1);
    for(int i = 0; i < m; i++)
    {
        int x, y;
        fin >> x >> y;
        adj[x].push_back(y);
        adj[y].push_back(x);
    }
    Graph g(7, adj);
    std::cout << g.dijkstra(1, 4);
    return 0;
}



// int Graph::dijkstra(int startNode, int endNode) {
//     std::vector<int> dist(numNodes, 1000000000);
//     std::priority_queue<std::pair<int, int>> pq;
//     dist[startNode] = 0;
//     pq.push(std::make_pair(0, startNode));
//     while(!pq.empty()) {
//         int currentNode = pq.top().second;
//         pq.pop();
//         for (int i = 0; i < adj[currentNode].size(); i++) {
//             int nextNode = adj[currentNode][i];
//             int weight = 1;
//             if (dist[nextNode] > dist[currentNode] + weight) {
//                 dist[nextNode] = dist[currentNode] + weight;
//                 pq.push(std::make_pair(-dist[nextNode], nextNode));
//             }
//         }
//     }
//     return dist[endNode];
// }

// Dijkstra's algorithm returning the shortest distance between two nodes
// int Graph::dijkstra(int startNode, int endNode) {
//     std::vector<edge> edges;
//     for (int i = 0; i < numNodes; i++) {
//         for (int j = 0; j < adj[i].size(); j++) {
//             edges.push_back({i, adj[i][j], 1});
//         }
//     }
//     return dijkstra(startNode, endNode, edges);
// }
