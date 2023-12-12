#include <iostream> 
#include <limits> 
#include <set>
#include <vector>
#include <iomanip> 
#include <functional>
#include <algorithm>
#include <random>
#include "utils.h" 

using namespace std;

// Utility classes and struct

struct VectorComparator {
    bool operator()(const vector<int>& a, const vector<int>& b) const {
        return a < b; // Use lexicographical comparison for vectors
    }
};

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

vector<int> generateRandomTour(int numCities){
    vector<int> tour(numCities);
    for (int i = 0; i < numCities; ++i) {
        tour[i] = i;
    }
    // Use a random device as a seed for the random number generator
    random_device rd;
    // Use the random device to seed the random number generator
    mt19937 gen(rd());
    // Shuffle the vector to generate a random permutation
    shuffle(tour.begin(), tour.end(), gen);
    return tour;
}

vector<vector<int>> initializePopulation(int N, int numCities){
    set<vector<int>, VectorComparator> generatedPopulationSet;
    for (int i = 0; generatedPopulationSet.size() < N; ++i) {
        vector<int> vec = generateRandomTour(numCities);
        generatedPopulationSet.insert(vec);
    }
    // Generate vector from set
    vector<vector<int>> generatedPopulation(generatedPopulationSet.begin(), generatedPopulationSet.end());
    return generatedPopulation;
}

void geneticAlgorithm(vector<vector<int>>& graph, vector<vector<int>>& edges, int startCity, int numCities, int numEdges){
    int populationSize;
    vector<int> initialSolution;
    vector<vector<int>> population, generatedRandomizedPopulation;

    populationSize = 25;
    // Generate inital solution using Krustal's MST algo
    initialSolution = kruskalMST(graph, edges, startCity, numCities, numEdges);

    // Initalize population
    population.push_back(initialSolution);
    generatedRandomizedPopulation = initializePopulation(populationSize-1,numCities);
    population.insert(population.end(), generatedRandomizedPopulation.begin(), generatedRandomizedPopulation.end());

    // Print generated populations
    std::cout << "Populations:" << std::endl;
    for (const auto& individual : population) {
        for (int value : individual) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    // Run generation loop:
    //      In every generation 
    //      loop for populationSize times: 
    //          perform selection, mutation and crossover -> generate offspring
    //      now we have 2*populationSize (half old and half new)
    //      choose populationSize best individuals from population using tour distance as metric and repeat above loop
    // Now find the best population amongnst current population
}





















