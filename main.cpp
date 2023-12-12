#include <iostream>
#include <fstream> 
#include <vector>
#include <limits> 
#include <cassert>
#include <iomanip> 
#include "utils.h"
using namespace std;

#define INF numeric_limits<int>::max()
 
int main()
{   
    int numCities, numEdges, startCity;

    // Read i/p from file
    ifstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        cerr << "Error opening input.txt" << endl;
        return 1;
    }
    inputFile >> numCities >> numEdges;

    assert(numEdges <= numCities*(numCities-1)/2);

    // Create a graph represented by an adjacency matrix
    vector<vector<int>> graph(numCities, vector<int>(numCities, INF));
    vector<vector<int>> edges;

    // Initialize diagonal elements to 0 (distance from a city to itself is 0)
    for (int i = 0; i < numCities; ++i) {
        graph[i][i] = 0;
    }

    // Read and process each edge
    for (int i = 0; i < numEdges; ++i) {
        int start, end, distance;
        inputFile >> start >> end >> distance;
        // add to graph
        graph[start][end] = distance;
        graph[end][start] = distance;
        // add to list of edges
        edges.push_back({start, end, distance});
        edges.push_back({end, start, distance});
    }

    // Close the file
    inputFile.close();

    // Print the adjacency matrix
    cout << "Adjacency Matrix:" << endl;
    printGraph(graph);

    // Let start_city = 0 by default
    startCity = 0;

    // Genetic algorithm
    geneticAlgorithm(graph, edges, startCity, numCities, numEdges);

    return 0;
}
