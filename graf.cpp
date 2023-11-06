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

// overloading the < operator for the edges
bool operator<(const edge &a, const edge &b)
{
    return a.weight < b.weight;
}

class Graph
{
private:
    int numNodes;
    std::vector<std::vector<int>> adj; // adjaency list
    std::vector<bool> visited;
    std::vector<edge> edges;
    std::vector<int> nodeValues; // node Values

    // functions for union find algorithm
    void makeSet(int node, std::vector<int> &parent, std::vector<int> &rank);
    int find(int node, std::vector<int> &parent);
    void unionSet(int node1, int node2, std::vector<int> &parent, std::vector<int> &rank);

    // function used by dfsCriticalConnections because I reset the visited vector first
    void dfsCritical(int node, int parentNode, std::vector<int> &lvl, std::vector<int> &low, std::vector<bool> &visited, std::vector<bool> &isArticulation, std::vector<std::vector<int>> &bridges);

    // returning the position of a node from a matrix indices(i and j)
    int getPos(int n, int x, int y) { return x * n + y; }

    // resetting the visited vector
    void resetVisited() { visited = std::vector<bool>(numNodes, false); }

public:
    Graph(int numNodes = 0);

    // initialising graph for shortest bridge problem, each element in the matrix is a node
    Graph(std::vector<std::vector<int>> matrix);

    Graph(int numNodes, std::vector<std::vector<int>> adj);
    Graph(int numNodes, std::vector<std::vector<int>> adj, std::vector<edge> edges);

    std::vector<std::vector<int>> getAdj() { return adj; }
    std::vector<int> getNodeValues() { return nodeValues; }
    std::vector<bool> getVisited() { return visited; }

    // usually most problems gives you a vector of edges; this method converts the vector into an adjaency list
    void constructConnections(std::vector<std::vector<int>> connections, bool isDirected = false);

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
    std::vector<edge> minimumSpanningTree();
};

Graph::Graph(int numNodes)
{
    this->numNodes = numNodes;
    adj.resize(numNodes);
    visited.resize(numNodes);
}

Graph::Graph(std::vector<std::vector<int>> matrix)
{
    numNodes = matrix.size() * matrix.size();
    adj.resize(numNodes);
    nodeValues.resize(numNodes);
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++)
        {
            if (matrix[i][j] == 1)
            {
                nodeValues[getPos(matrix.size(), i, j)] = 1;
            }
            if (i - 1 >= 0)
                adj[getPos(matrix.size(), i, j)].push_back((i - 1) * matrix.size() + j);
            if (i + 1 < matrix.size())
                adj[getPos(matrix.size(), i, j)].push_back((i + 1) * matrix.size() + j);
            if (j - 1 >= 0)
                adj[getPos(matrix.size(), i, j)].push_back(getPos(matrix.size(), i, j) - 1);
            if (j + 1 < matrix.size())
                adj[getPos(matrix.size(), i, j)].push_back(getPos(matrix.size(), i, j) + 1);
        }
    }
}

// initialising graph using the number of nodes and the adjacency list
Graph::Graph(int numNodes, std::vector<std::vector<int>> adj)
{
    this->numNodes = numNodes;
    this->adj = adj;
    this->visited = std::vector<bool>(numNodes, false);
}

Graph::Graph(int numNodes, std::vector<std::vector<int>> adj, std::vector<edge> edges)
{
    this->numNodes = numNodes;
    this->adj = adj;
    this->edges = edges;
    this->visited = std::vector<bool>(numNodes, false);
}

void Graph::constructConnections(std::vector<std::vector<int>> connections, bool isDirected)
{
    if (!isDirected)
    {
        for (auto edge : connections)
        {
            adj[edge[0]].push_back(edge[1]);
            adj[edge[1]].push_back(edge[0]);
            edges.push_back({edge[0], edge[1], 1});
        }
        return;
    }
    for (auto edge : connections)
    {
        adj[edge[0]].push_back(edge[1]);
        edges.push_back({edge[0], edge[1], 1});
    }
}

void Graph::addNode()
{
    numNodes++;
}

void Graph::addEdge(int node1, int node2, int weight, bool isDirected)
{
    if (node1 == node2)
        return;

    // Check if the edge already exists in the adjacency list for undirected graph
    if (!isDirected)
    {
        for (const auto &neighbor : adj[node1])
        {
            if (neighbor == node2)
            {
                // Edge already exists, so skip adding it again
                return;
            }
        }
    }
    // Add the edge to the adjacency list
    adj[node1].push_back(node2);
    edges.push_back({node1, node2, weight});

    // If the graph is undirected, add the reverse edge
    if (!isDirected)
    {
        adj[node2].push_back(node1);
    }
}

std::vector<int> Graph::dfs(int startNode)
{
    resetVisited();
    std::stack<int> s;
    std::vector<int> returningDFS; // the vector of the DFS that will be returned
    visited[startNode] = true;
    s.push(startNode);
    while (!s.empty())
    {
        int currentNode = s.top();
        s.pop();
        returningDFS.push_back(currentNode); // adding the current node to the DFS
        for (int i = 0; i < adj[currentNode].size(); i++)
        {
            int neighbour = adj[currentNode][i];
            if (!visited[neighbour])
            {
                visited[neighbour] = true;
                s.push(neighbour);
            }
        }
    }
    return returningDFS;
}

std::pair<std::vector<bool>, std::vector<std::vector<int>>> Graph::dfsCriticalConnections(int node)
{
    resetVisited();
    int parentNode = -1;
    std::vector<int> lvl(numNodes, -1);
    std::vector<int> low(numNodes, -1);
    std::vector<bool> isArticulation(numNodes, false);
    std::vector<std::vector<int>> bridges(numNodes);
    dfsCritical(node, parentNode, lvl, low, visited, isArticulation, bridges);
    for (int i = 0; i < numNodes; ++i)
    {
        if (!visited[i])
        {
            dfsCritical(i, -1, lvl, low, visited, isArticulation, bridges);
        }
    }
    return {isArticulation, bridges};
    // std::cout << "Articulation points:\n";
    // for(int i = 0; i < numNodes; ++i) {
    //     if(isArticulation[i]) {
    //         std::cout << i << " ";
    //     }
    // }
    // std::cout << "\n Bridges:\n";
    // for(int i = 0; i < bridges.size(); ++i) {
    //     for(int j = 0; j < bridges[i].size(); ++j) {
    //         std::cout << bridges[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }
}

/**
 * @brief Perform Depth-First Search to identify critical nodes (articulation points) and bridges
 * in an undirected graph.
 *
 * @param node Current node being explored.
 * @param parentNode Parent node of the current node in the DFS traversal.
 * @param lvl Vector to store levels of nodes in DFS tree.
 * @param low Vector to store lowest level reachable from the current node.
 * @param visited Vector to mark nodes as visited during DFS.
 * @param isArticulation Vector to mark nodes as articulation points.
 * @param bridges Vector of vectors to store bridge edges (articulation points).
 */
void Graph::dfsCritical(int node, int parentNode, std::vector<int> &lvl, std::vector<int> &low, std::vector<bool> &visited, std::vector<bool> &isArticulation, std::vector<std::vector<int>> &bridges)
{
    // Set level and low value for the current node
    if (parentNode == -1)
        lvl[node] = 1;
    else
        lvl[node] = lvl[parentNode] + 1;
    low[node] = lvl[node];

    visited[node] = true;

    // Count children of the current node in the DFS tree
    int children = 0;

    // Explore neighbors of the current node
    for (auto neighbour : adj[node])
    {
        if (!visited[neighbour])
        {
            // Recur for unvisited neighbors
            dfsCritical(neighbour, node, lvl, low, visited, isArticulation, bridges);
            children++;

            // Update low value of the current node based on its children
            low[node] = std::min(low[node], low[neighbour]);

            // Check for articulation points (nodes with no back edge)
            if (low[neighbour] >= lvl[node] && node != 0)
            {
                isArticulation[node] = true;
            }

            // Check for bridges (edges connecting nodes with no back edge)
            if (low[neighbour] > lvl[node])
            {
                bridges.push_back({node, neighbour});
            }
        }
        else if (neighbour != parentNode)
        {
            // Update low value if the neighbour is not the parent (back edge)
            low[node] = std::min(low[node], lvl[neighbour]);
        }
    }

    // Check for root node with multiple children
    if (parentNode == -1 && children > 1)
        isArticulation[node] = true;
}

std::vector<int> Graph::bfs(int startNode)
{
    resetVisited();

    // the vector of the BFS that will be returned
    std::vector<int> returningBFS;

    std::queue<int> q;
    visited[startNode] = true;
    q.push(startNode);

    while (!q.empty())
    {
        int currentNode = q.front();
        q.pop();
        returningBFS.push_back(currentNode); // adding the current node to the BFS
        for (int i = 0; i < adj[currentNode].size(); i++)
        {
            if (!visited[adj[currentNode][i]])
            {
                visited[adj[currentNode][i]] = true;
                q.push(adj[currentNode][i]);
            }
        }
    }
    return returningBFS;
}

int Graph::dijkstra(std::vector<int> startNodes, std::vector<int> endNodes)
{
    std::vector<int> dist(numNodes, 1000000000); // vector of distances from the start node to the other nodes
    std::priority_queue<std::pair<int, int>> pq; // priority queue of pairs of distances and nodes
    if (edges.empty())                           // if the vector of edges is empty then we create it
        for (int i = 0; i < adj.size(); ++i)
        {
            for (int j = 0; j < adj[i].size(); ++j)
            {
                edges.push_back({i, adj[i][j], 1});
            }
        }
    for (auto startNode : startNodes)
    {
        dist[startNode] = 0;
        pq.push(std::make_pair(0, startNode));
    }
    while (!pq.empty())
    {
        int currentNode = pq.top().second;
        pq.pop();
        for (int i = 0; i < edges.size(); i++)
        {
            // finding the nodes adjacents to the current node
            if (edges[i].node1 == currentNode)
            {
                int nextNode = edges[i].node2;
                int weight = edges[i].weight;

                // if the distance from the next node is bigger than the distance from the current node + the weight of the edge
                if (dist[nextNode] > dist[currentNode] + weight)
                {
                    // the distance from the next node becomes the distance from the current node + the weight of the edge
                    dist[nextNode] = dist[currentNode] + weight;
                    pq.push(std::make_pair(-dist[nextNode], nextNode));
                }
            }
        }
    }
    int minDist = 1000000000;

    // search for the shortest path between the start nodes and the end nodes
    for (auto endNode : endNodes)
    {
        if (dist[endNode] < minDist)
            minDist = dist[endNode];
    }
    return minDist; // returning the distance from the start node to the end node
}

std::vector<int> Graph::bellmanFord(int startNode)
{
    std::vector<int> dist(numNodes, 1000000000);
    dist[startNode] = 0;

    // for every node except the start node
    for (int i = 0; i < numNodes - 1; i++)
    {
        for (int j = 0; j < edges.size(); j++)
        {
            int node1 = edges[j].node1;
            int node2 = edges[j].node2;
            int weight = edges[j].weight;
            if (dist[node2] > dist[node1] + weight)
            {                                       // if the distance from the second node is bigger than the distance from the first node + the weight of the edge then:
                dist[node2] = dist[node1] + weight; // the distance from the second node becomes the distance from the first node + the weight of the edge
            }
        }
    }
    return dist; // returning the vector of distances
}

// check if the graph is bipartite
bool Graph::isBipartite()
{
    int node = 0;
    std::queue<int> q;
    q.push(node);
    std::vector<int> colors(numNodes, -1); // vector of colors for the nodes (-1 = not colored, 0 = color 1, 1 = color 2)

    // check for every component of the graph
    for (int i = 0; i < numNodes; i++)
    {
        // if the node is not colored then:
        if (colors[i] == -1)
        {
            colors[i] = 0; // color the node with color 1
            q.push(i);
        }

        // while the queue is not empty means that we still have nodes to color
        while (!q.empty())
        {
            int currentNode = q.front();
            q.pop();
            for (int i = 0; i < adj[currentNode].size(); ++i)
            {
                int nextNode = adj[currentNode][i];

                if (colors[nextNode] == -1)
                {
                    // color the next node with the opposite color of the current node
                    colors[nextNode] = 1 - colors[currentNode];
                    q.push(nextNode);
                }

                // if the next node is colored with the same color as the current node then
                else if (colors[nextNode] == colors[currentNode])
                {
                    // the graph is not bipartite
                    return false;
                }
            }
        }
    }

    // if the queue is empty means that we colored all the nodes and the graph is bipartite
    return true;
}

std::vector<int> Graph::topologicalSort()
{
    std::vector<int> returningTopologicalSort; // the vector of the topological sort that will be returned
    std::vector<int> inDegree(numNodes, 0);
    for (int i = 0; i < numNodes; ++i)
    {
        for (int j = 0; j < adj[i].size(); ++j)
        {
            inDegree[adj[i][j]]++;
        }
    }

    // queue of nodes with in degree 0
    std::queue<int> q;

    for (int i = 0; i < numNodes; ++i)
    {
        if (inDegree[i] == 0)
        {
            q.push(i);
        }
    }

    while (!q.empty())
    {
        int currentNode = q.front();
        q.pop();
        returningTopologicalSort.push_back(currentNode);

        for (int i = 0; i < adj[currentNode].size(); ++i)
        {
            int nextNode = adj[currentNode][i];

            // decreasing the in degree of the next node because we removed the current node from the graph
            inDegree[nextNode]--;
            if (inDegree[nextNode] == 0)
            {
                q.push(nextNode);
            }
        }
    }
    return returningTopologicalSort;
}

// UNION FIND ALGORITHM
void Graph::makeSet(int node, std::vector<int> &parent, std::vector<int> &rank)
{
    // initialising every node as a parent of itself
    parent[node] = node;
    rank[node] = 0;
}

int Graph::find(int node, std::vector<int> &parent)
{
    // going up in the tree until the node is the parent of itself
    while (node != parent[node])
        node = parent[node];
    return node;
}

void Graph::unionSet(int node1, int node2, std::vector<int> &parent, std::vector<int> &rank)
{
    // finding the roots
    int root1 = find(node1, parent);
    int root2 = find(node2, parent);

    // if the roots are the same then the nodes are in the same set
    if (root1 == root2)
        return;

    if (rank[root1] > rank[root2])
        parent[root2] = root1;
    else
    {
        parent[root1] = root2;
        if (rank[root1] == rank[root2])
            rank[root2]++;
    }
}

std::vector<edge> Graph::minimumSpanningTree()
{
    std::vector<edge> returningMST; // the vector of edges that will be returned
    std::vector<edge> sortedEdges = edges;
    std::sort(sortedEdges.begin(), sortedEdges.end()); // sorting the edges by weight
    std::vector<int> parent(numNodes);
    std::vector<int> rank(numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        // initialising the parents and ranks of the nodes
        makeSet(i, parent, rank);
    }
    for (auto edge : sortedEdges)
    {
        int node1 = edge.node1;
        int node2 = edge.node2;

        // if the nodes are not in the same set then:
        if (find(node1, parent) != find(node2, parent))
        {
            // add the edge to the minimum spanning tree
            returningMST.push_back(edge);

            // merge the sets
            unionSet(node1, node2, parent, rank);
        }
    }
    return returningMST;
}

int main()
{
    std::ifstream fin("graf.in");
    // std::ifstream fin("graf2.in");
    int n, m;
    fin >> n >> m;
    std::vector<std::vector<int>> adj(n + 1);
    for (int i = 0; i < m; i++)
    {
        int x, y;
        fin >> x >> y;
        adj[x].push_back(y);
        adj[y].push_back(x);
    }
    Graph g(n, adj);
    g.dfsCriticalConnections(0);
    return 0;
}

// int main () {
//     std::ifstream fin("graf.in");
//     // std::ifstream fin("graf2.in");
//     int n, m;
//     fin >> n >> m;
//     // std::vector<std::vector<int>> adj(n+1);
//     // for(int i = 0; i < m; i++)
//     // {
//     //     int x, y;
//     //     fin >> x >> y;
//     //     adj[x].push_back(y);
//     //     adj[y].push_back(x);
//     // }
//     std::vector<std::vector<int>> adj(n);
//     std::vector<edge> edges(m);
//     for(int i = 0; i < m; i++)
//     {
//         int x, y, w;
//         fin >> x >> y >> w;
//         adj[x].push_back(y);
//         adj[y].push_back(x);
//         edges[i] = {x, y, w};
//     }
//     Graph g(n, adj, edges);
//     // std::cout << g.dijkstra(1, 4);
//     std::vector<edge> v = g.minimumSpanningTree();
//     for(int i = 0; i < v.size(); ++i) {
//         std::cout << v[i].node1 << " " << v[i].node2 << " " << v[i].weight << "\n";
//     }
//     return 0;
// }