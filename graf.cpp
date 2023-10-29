#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <fstream>
#include <unordered_set>

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
    std::vector<edge> bridges;

    // reseting class vectors
    void resetVisited();
    void resetLvl();
    void resetLow();
    void resetIsArticulation();
    void resetAll();

    //getters
    const std::vector<bool> getVisited() const { return visited; }
    const std::vector<bool> getArticulation() const { return isArticulation; }
    const std::vector<edge> getBridge() const { return bridges; }
    
    // functions for union find algorithm
    void makeSet(int node, std::vector<int>& parent, std::vector<int>& rank);
    int find(int node, std::vector<int>& parent);
    void unionSet(int node1, int node2, std::vector<int>& parent, std::vector<int>& rank);

public:
    Graph(int numNodes = 0);
    Graph(int numNodes, std::vector<std::vector<int>> adj);
    Graph(int numNodes, std::vector<std::vector<int>> adj, std::vector<edge> edges);

    void addEdge(int node1, int node2, int weight = 1, bool isDirected = false);
    void addNode();
    std::vector<int> dfs(int startNode);
    void dfsCriticalConditions(int node, int parentNode = 0);
    std::vector<int> bfs(int startNode);
    int dijkstra(int startNode, int endNode);
    std::vector<int> bellmanFord(int startNode);
    bool isBipartite(int node, int parentNode);
    std::vector<int> topologicalSort();
    std::vector<edge> minimumSpanningTree();
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
    resetAll(); // reseting the visited vector
    std::stack<int> s; // stack of nodes
    std::vector<int> returningDFS; // the vector of the DFS that will be returned
    visited[startNode] = true; // marking the start node as visited
    s.push(startNode); // adding the start node to the stack
    returningDFS.push_back(startNode); // adding the start node to the DFS
    while(!s.empty()) { // while the stack is not empty means that we still have nodes to visit
        int currentNode = s.top(); // the current node is the top node in the stack
        s.pop(); // removing the current node from the stack
        returningDFS.push_back(currentNode); // adding the current node to the DFS
        for (int i = 0; i < adj[currentNode].size(); i++) {
            if (!visited[adj[currentNode][i]]) { // if the node is not visited then:
                visited[adj[currentNode][i]] = true; // mark the node as visited
                s.push(adj[currentNode][i]); // add the node to the stack
            }
        }
    }
    return returningDFS;
}

// DFS implementation for finding articulation points and bridges -> Articulation points will be in isArticulation from class
// and bridges will be in the returningBridges vector
void Graph::dfsCriticalConditions(int node, int parentNode) {
    resetAll(); // reseting the visited vector
    lvl[node] = lvl[parentNode]  + 1;
    low[node] = lvl[node];
    visited[node] = 1;
    int children = 0;
    for(auto neighbour: adj[node]) {
        if(!visited[neighbour]) {
            dfsCriticalConditions(neighbour, node);
            children++;
            low[node] = std::min(low[node], low[neighbour]);
            if(low[neighbour] >= lvl[node] && node != 0) {
                isArticulation[node] = true;
            }
            if(low[neighbour] > lvl[node]) {
                bridges.push_back({node, neighbour, 1});
            }
        } else if(neighbour != parentNode) {
            low[node] = std::min(low[node], lvl[neighbour]);
        }
    }
    if(parentNode == 0 && children > 1)
        isArticulation[node] = true;
}

// Breadth first search returning a vector with the order of the nodes visited
std::vector<int> Graph::bfs(int startNode) {
    resetAll(); // reseting the visited vector
    std::vector<int> returningBFS; // the vector of the BFS that will be returned
    std::queue<int> q; // queue of nodes
    visited[startNode] = true; // marking the start node as visited
    q.push(startNode); // adding the start node to the queue
    while(!q.empty()) { // while the queue is not empty means that we still have nodes to visit
        int currentNode = q.front(); // the current node is the first node in the queue
        q.pop(); // removing the current node from the queue
        returningBFS.push_back(currentNode); // adding the current node to the BFS
        for (int i = 0; i < adj[currentNode].size(); i++) { // for every node adjacent to the current node
            if (!visited[adj[currentNode][i]]) { // if the node is not visited then:
                visited[adj[currentNode][i]] = true; // mark the node as visited
                q.push(adj[currentNode][i]); // add the node to the queue
            }
        }
    }
    return returningBFS;
}

// TO DO: Dijkstra algorithm should return VECTOR not INT and you can get the SIZE by VECTOR.SIZE() instead
// Dijkstra's algorithm returning the shortest distance between two nodes using a vector of edges
int Graph::dijkstra(int startNode, int endNode) {
    std::vector<int> dist(numNodes, 1000000000); // vector of distances from the start node to the other nodes
    std::priority_queue<std::pair<int, int>> pq; // priority queue of pairs of distances and nodes
    if(edges.empty()) // if the vector of edges is empty then we create it
        for(int i = 0; i < adj.size(); ++i) { 
            for(int j = 0; j < adj[i].size(); ++j) {
                edges.push_back({i, adj[i][j], 1}); // adding the edges to the vector
            }
    }
    dist[startNode] = 0; // the distance from the start node to itself is 0
    pq.push(std::make_pair(0, startNode)); // adding the start node to the priority queue
    while(!pq.empty()) {
        int currentNode = pq.top().second; // the current node is the second element of the pair
        pq.pop(); // removing the current node from the priority queue
        for (int i = 0; i < edges.size(); i++) { // for every edge
            if (edges[i].node1 == currentNode) { // if the first node of the edge is the current node then:
                int nextNode = edges[i].node2; // the next node is the second node of the edge
                int weight = edges[i].weight;
                if (dist[nextNode] > dist[currentNode] + weight) { // if the distance from the next node is bigger than the distance from the current node + the weight of the edge then:
                    dist[nextNode] = dist[currentNode] + weight; // the distance from the next node becomes the distance from the current node + the weight of the edge
                    pq.push(std::make_pair(-dist[nextNode], nextNode)); // adding the next node to the priority queue
                }
            }
        }
    }
    return dist[endNode]; // returning the distance from the start node to the end node
}

std::vector<int> Graph::bellmanFord(int startNode) {
    std::vector<int> dist(numNodes, 1000000000); // vector of distances from the start node to the other nodes
    dist[startNode] = 0; // the distance from the start node to itself is 0
    for (int i = 0; i < numNodes - 1; i++) { // for every node except the start node
        for (int j = 0; j < edges.size(); j++) { // for every edge
            int node1 = edges[j].node1; // first node of the edge
            int node2 = edges[j].node2; // second node of the edge
            int weight = edges[j].weight; // weight of the edge
            if (dist[node2] > dist[node1] + weight) { // if the distance from the second node is bigger than the distance from the first node + the weight of the edge then:
                dist[node2] = dist[node1] + weight; // the distance from the second node becomes the distance from the first node + the weight of the edge
            }
        }
    }
    return dist; // returning the vector of distances
}

// check if the graph is bipartite, before using it you should call ResetAll() method because the visited vector might be changed
bool Graph::isBipartite(int node, int parentNode) {
    std::queue<int> q; // queue of nodes
    q.push(node); // adding the first node to the queue
    std::vector<int> colors(numNodes, -1); // vector of colors for the nodes (-1 = not colored, 0 = color 1, 1 = color 2)
    colors[node] = 0; // coloring the first node with color 1
    while(!q.empty()) { // while the queue is not empty means that we still have nodes to color
        int currentNode = q.front();
        q.pop();
        for(int i = 0; i < adj[currentNode].size(); ++i) {
            int nextNode = adj[currentNode][i]; // the next node is the node that we are going to color
            if(colors[nextNode] == -1) { // if the next node is not colored then:
                colors[nextNode] = 1 - colors[currentNode]; // color the next node with the opposite color of the current node
                q.push(nextNode); // add the next node to the queue
            } else if(colors[nextNode] == colors[currentNode]) { // if the next node is colored with the same color as the current node then:
                return false; // the graph is not bipartite
            }
        }
    }
    return true; // if the queue is empty means that we colored all the nodes and the graph is bipartite
}

// topological sorting of a graph using Kahn's algorithm returning a vector with the order of the nodes
std::vector<int> Graph::topologicalSort() {
    std::vector<int> returningTopologicalSort; // the vector of the topological sort that will be returned
    std::vector<int> inDegree(numNodes, 0); // the vector of the in degrees of the nodes
    for(int i = 0; i < numNodes; ++i) {
        for(int j = 0; j < adj[i].size(); ++j) {
            inDegree[adj[i][j]]++; // calculating the in degrees of the nodes
        }
    }
    std::queue<int> q; // queue of nodes with in degree 0
    for(int i = 0; i < numNodes; ++i) {
        if(inDegree[i] == 0) {
            q.push(i); // adding the nodes with in degree 0 to the queue
        }
    }
    while(!q.empty()) { // while the queue is not empty means that we still have nodes to add to the topological sort
        int currentNode = q.front(); // the current node is the first node in the queue
        q.pop(); // removing the current node from the queue
        returningTopologicalSort.push_back(currentNode); // adding the current node to the topological sort
        for(int i = 0; i < adj[currentNode].size(); ++i) {
            int nextNode = adj[currentNode][i]; 
            inDegree[nextNode]--; // decreasing the in degree of the next node because we removed the current node from the graph
            if(inDegree[nextNode] == 0) { 
                q.push(nextNode); // if the in degree of the next node is 0 then we add it to the queue
            }
        }
    }
    return returningTopologicalSort;
}

// UNION FIND ALGORITHM
void Graph::makeSet(int node, std::vector<int>& parent, std::vector<int>& rank)
{
    parent[node] = node; // initialising every node as a parent of itself
    rank[node] = 0; // initialising every node with rank 0
}

int Graph::find(int node, std::vector<int>& parent)
{
    while(node != parent[node]) // finding the root of the node
        node = parent[node]; // going up in the tree until the node is the parent of itself
    return node; // returning the root of the node
}

void Graph::unionSet(int node1, int node2, std::vector<int>& parent, std::vector<int>& rank)
{
    int root1 = find(node1, parent); // finding the root of the first node
    int root2 = find(node2, parent); // finding the root of the second node
    if(root1 == root2) // if the roots are the same then the nodes are in the same set
        return;
    if(rank[root1] > rank[root2]) // if the rank of the first root is bigger than the rank of the second root then: 
        parent[root2] = root1; // the first root becomes the parent of the second root
    else {
        parent[root1] = root2; 
        if(rank[root1] == rank[root2]) // if the ranks are equal then the rank of the second root increases by 1
            rank[root2]++;
    }
}


// minimum spanning tree using Kruskal's algorithm returning a vector of edges
std::vector<edge> Graph::minimumSpanningTree() {
    std::vector<edge> returningMST; // the vector of edges that will be returned
    std::vector<edge> sortedEdges = edges;
    std::sort(sortedEdges.begin(), sortedEdges.end()); // sorting the edges by weight
    std::vector<int> parent(numNodes); // vector of parents for the union find algorithm
    std::vector<int> rank(numNodes); // vector of ranks for the union find algorithm
    for(int i = 0; i < numNodes; i++) {
        makeSet(i, parent, rank); // initialising the parents and ranks of the nodes
    }
    for(auto edge : sortedEdges) {
        int node1 = edge.node1; // first node of the edge
        int node2 = edge.node2; // second node of the edge
        if(find(node1, parent) != find(node2, parent)) // if the nodes are not in the same set then:
        {
            returningMST.push_back(edge); // add the edge to the minimum spanning tree
            unionSet(node1, node2, parent, rank); // union the sets of the nodes
        }
    }
    return returningMST;
}

// articulation points and bridges


// deseneaza pe o foaie Graphuri ca sa intelegi articulation points and bridges
    
int main () {
    std::ifstream fin("graf.in");
    // std::ifstream fin("graf2.in");
    int n, m;
    fin >> n >> m;
    // std::vector<std::vector<int>> adj(n+1);
    // for(int i = 0; i < m; i++)
    // {
    //     int x, y;
    //     fin >> x >> y;
    //     adj[x].push_back(y);
    //     adj[y].push_back(x);
    // }
    std::vector<std::vector<int>> adj(n);
    std::vector<edge> edges(m);
    for(int i = 0; i < m; i++)
    {
        int x, y, w;
        fin >> x >> y >> w;
        adj[x].push_back(y);
        adj[y].push_back(x);
        edges[i] = {x, y, w};
    }
    Graph g(n, adj, edges);
    // std::cout << g.dijkstra(1, 4);
    std::vector<edge> v = g.minimumSpanningTree();
    for(int i = 0; i < v.size(); ++i) {
        std::cout << v[i].node1 << " " << v[i].node2 << " " << v[i].weight << "\n";
    }
    return 0;
}