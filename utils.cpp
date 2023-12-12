#include <iostream> 
#include <limits> 
#include <vector>
#include <iomanip> 
#include <functional>
#include "utils.h" 

using namespace std;

// Utility classes

class DisjointSet {
private:
    vector<int> parent, rank;

public:
    DisjointSet(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    // Find the set (representative) of the given element
    int findSet(int x) {
        if (x != parent[x]) {
            parent[x] = findSet(parent[x]);
        }
        return parent[x];
    }

    // Merge two sets if they are different
    void unionSets(int x, int y) {
        int rootX = findSet(x);
        int rootY = findSet(y);

        if (rootX != rootY) {
            if (rank[rootX] < rank[rootY]) {
                swap(rootX, rootY);
            }
            parent[rootY] = rootX;
            if (rank[rootX] == rank[rootY]) {
                rank[rootX]++;
            }
        }
    }
};

// Utility function definitions

void printGraph(vector<vector<int>>& graph){
    for (auto& row : graph) {
        for (int distance : row) {
            if (distance == INF) {
                cout << setw(5) << "INF";
            } else {
                cout << setw(5) << distance;
            }
        }
        cout << endl;
    }
}

vector<int> kruskalMST(vector<vector<int>>& graph, vector<vector<int>>& edges, int startCity, int numCities, int numEdges) {
    int n = numCities;

    // Create disjoint set for each city
    DisjointSet ds(n);

    // Kruskal's algorithm to find MST
    vector<vector<int>> mstEdges;
    for (auto& edge : edges) {
        int u = edge[0];
        int v = edge[1];

        // Check if adding this edge forms a cycle in the MST
        if (ds.findSet(u) != ds.findSet(v)) {
            // Add the edge to the MST
            mstEdges.push_back(edge);
            // Merge the sets of vertices u and v
            ds.unionSets(u, v);
        }
    }

    // Create an adjacency list to represent the MST
    vector<vector<int>> adjacencyList(n);
    for (auto& edge : mstEdges) {
        adjacencyList[edge[0]].push_back(edge[1]);
        adjacencyList[edge[1]].push_back(edge[0]);
    }

    // Perform a depth-first traversal to create a TSP path
    vector<int> tspPath;
    vector<bool> visited(n, false);

    // Depth-first traversal function
    function<void(int)> dfs = [&](int current) {
        tspPath.push_back(current);
        visited[current] = true;
        for (int neighbor : adjacencyList[current]) {
            if (!visited[neighbor]) {
                dfs(neighbor);
            }
        }
    };
    dfs(startCity);
    return tspPath;
}

void geneticAlgorithm(vector<vector<int>>& graph, vector<vector<int>>& edges, int startCity, int numCities, int numEdges){
    vector<int> initialSolution;
    initialSolution = kruskalMST(graph, edges, startCity, numCities, numEdges);
    for (size_t i = 0; i < initialSolution.size(); ++i) {
        cout << initialSolution[i] << " ";
    }
    cout << endl;
}





















