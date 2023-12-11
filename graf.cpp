#include "graf.h"

// overloading the < operator for the edges
bool operator<(const edge &a, const edge &b)
{
    return a.weight < b.weight;
}

Graph::Graph(int numNodes)
{
    this->numNodes = numNodes;
    adj.resize(numNodes);
    visited.resize(numNodes);
}

// initialising graph using the number of nodes and the adjacency list
Graph::Graph(std::vector<std::vector<std::pair<int, int>>> adj)
{
    this->numNodes = adj.size();
    this->adj = adj;
    this->visited = std::vector<bool>(numNodes, false);
}

Graph::Graph(int numNodes, std::vector<std::vector<int>> connections, bool isDirected)
{
    this->numNodes = numNodes;
    if (!isDirected)
    {
        for (auto edge : connections)
        {

            // if the connections vector has edge1, edge2 and the cost between them, push the cost as well
            if (edge.size() > 2)
            {
                adj[edge[0]].push_back({edge[1], edge[2]});
                adj[edge[1]].push_back({edge[0], edge[2]});
            }
            // else push 1 for the cost
            else
            {
                adj[edge[0]].push_back({edge[1], 1});
                adj[edge[1]].push_back({edge[0], 1});
            }
        }
        return;
    }
    for (auto edge : connections)
    {
        if (edge.size() > 2)
            adj[edge[0]].push_back({edge[1], edge[2]});
        else
            adj[edge[0]].push_back({edge[1], 1});
    }
}

void Graph::addNode()
{
    numNodes++;
}

void Graph::addEdge(int node1, int node2, int weight, bool isDirected)
{
    if (node1 == node2 || node1 > numNodes || node2 > numNodes)
        return;

    // Check if the edge already exists in the adjacency list for undirected graph

    for (const auto &neighbor : adj[node1])
    {
        if (neighbor.first == node2)
        {
            // Edge already exists, so skip adding it again
            return;
        }
    }

    // Add the edge to the adjacency list
    adj[node1].push_back({node2, weight});

    // If the graph is undirected, add the reverse edge
    if (!isDirected)
    {
        adj[node2].push_back({node1, weight});
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
            int neighbour = adj[currentNode][i].first;
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
    std::vector<std::vector<int>> bridges;
    dfsCritical(node, parentNode, lvl, low, isArticulation, bridges);
    for (int i = 0; i < numNodes; ++i)
    {
        if (!visited[i])
        {
            dfsCritical(i, -1, lvl, low, isArticulation, bridges);
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
void Graph::dfsCritical(int node, int parentNode, std::vector<int> &lvl, std::vector<int> &low, std::vector<bool> &isArticulation, std::vector<std::vector<int>> &bridges)
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
        // neighbour is a pair {node, weight} => neigh is the node
        int neigh = neighbour.first;
        if (!visited[neigh])
        {
            // Recur for unvisited neighbors
            dfsCritical(neigh, node, lvl, low, isArticulation, bridges);
            children++;

            // Update low value of the current node based on its children
            low[node] = std::min(low[node], low[neigh]);

            // Check for articulation points (nodes with no back edge)
            if (low[neigh] >= lvl[node] && node != 0)
            {
                isArticulation[node] = true;
            }

            // Check for bridges (edges connecting nodes with no back edge)
            if (low[neigh] > lvl[node])
            {
                bridges.push_back({node, neigh});
            }
        }
        else if (neigh != parentNode)
        {
            // Update low value if the neighbour is not the parent (back edge)
            low[node] = std::min(low[node], lvl[neigh]);
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
            int neigh = adj[currentNode][i].first;
            if (!visited[neigh])
            {
                visited[neigh] = true;
                q.push(neigh);
            }
        }
    }
    return returningBFS;
}

int Graph::dijkstra(std::vector<int> startNodes, std::vector<int> endNodes)
{
    std::vector<int> dist(numNodes, 1000000000); // vector of distances from the start node to the other nodes
    std::priority_queue<std::pair<int, int>> pq; // priority queue of pairs of distances and nodes
    for (auto startNode : startNodes)
    {
        dist[startNode] = 0;
        pq.push(std::make_pair(0, startNode));
    }
    while (!pq.empty())
    {
        int currentNode = pq.top().second;
        pq.pop();
        for (auto neigh : adj[currentNode])
        {
            int nextNode = neigh.first;
            int weight = neigh.second;

            // if the distance from the next node is bigger than the distance from the current node + the weight of the edge
            if (dist[nextNode] > dist[currentNode] + weight)
            {
                // the distance from the next node becomes the distance from the current node + the weight of the edge
                dist[nextNode] = dist[currentNode] + weight;
                pq.push(std::make_pair(-dist[nextNode], nextNode));
            }
        }
    }
    int minDist = 1000000000;

    // search for the shortest path between the start nodes and the end nodes
    for (auto endNode : endNodes)
        minDist = std::min(minDist, dist[endNode]);

    return minDist; // returning the distance from the start node to the end node
}

std::vector<int> Graph::bellmanFord(int startNode)
{
    std::vector<int> dist(numNodes, INT_MAX);
    dist[startNode] = 0;

    std::vector<edge> edges;

    // inserting all edges into a single vector
    for (int i = 0; i < numNodes; ++i)
    {
        for (auto neigh : adj[i])
        {
            edges.push_back({i, neigh.first, neigh.second});
        }
    }

    // relaxing all edges numNodes - 1 times
    for (int i = 1; i <= numNodes - 1; i++)
    {
        for (int j = 0; j < edges.size(); j++)
        {
            int node1 = edges[j].node1;
            int node2 = edges[j].node2;
            int weight = edges[j].weight;
            if (dist[node1] != INT_MAX && dist[node2] > dist[node1] + weight)
            {                                       // if the distance from the second node is bigger than the distance from the first node + the weight of the edge then:
                dist[node2] = dist[node1] + weight; // the distance from the second node becomes the distance from the first node + the weight of the edge
            }
        }
    }

    // checking for negative-weight cycles
    for (int i = 0; i < edges.size(); i++)
    {
        int node1 = edges[i].node1;
        int node2 = edges[i].node2;
        int weight = edges[i].weight;
        if (dist[node1] != INT_MAX && dist[node2] > dist[node1] + weight)
        {
            std::cout << "Graph contains negative weight cycle";
            return {};
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
                int nextNode = adj[currentNode][i].first;

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
            inDegree[adj[i][j].first]++;
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
            int nextNode = adj[currentNode][i].first;

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

std::vector<edge> Graph::kruskalMST()
{
    std::vector<edge> returningMST; // the vector of edges that will be returned
    std::vector<edge> sortedEdges;
    for (int i = 0; i < numNodes; ++i)
        for (auto neigh : adj[i])
            sortedEdges.push_back({i, neigh.first, neigh.second});
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

std::vector<edge> Graph::primMST()
{
    int startNode = 0;
    std::vector<edge> returningMST; // the vector of edges that will be returned
    std::vector<int> dist(numNodes, 1e9);
    std::vector<int> parent(numNodes, -1);
    std::vector<bool> visited(numNodes, false);
    dist[startNode] = 0;
    for (int i = 0; i < numNodes; ++i)
    {
        int currentNode = -1;
        for (int j = 0; j < numNodes; ++j)
        {
            if (!visited[j] && (currentNode == -1 || dist[j] < dist[currentNode]))
            {
                currentNode = j;
            }
        }
        if (dist[currentNode] == 1e9)
        {
            std::cout << "No MST!";
            return {};
        }
        visited[currentNode] = true;
        if (parent[currentNode] != -1)
        {
            returningMST.push_back({parent[currentNode], currentNode, dist[currentNode]});
        }
        for (auto neigh : adj[currentNode])
        {
            int nextNode = neigh.first;
            int weight = neigh.second;
            if (dist[nextNode] > weight)
            {
                dist[nextNode] = weight;
                parent[nextNode] = currentNode;
            }
        }
    }
    return returningMST;
}