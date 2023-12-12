#include <iostream>
#include <fstream> 
#include <vector>
#include <limits> 
#include <iomanip> 
#include "utils.h"
using namespace std;

#define INF numeric_limits<int>::max()
 
int main()
{   
    int numCities, numEdges;

    // Read i/p from file
    ifstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        cerr << "Error opening input.txt" << endl;
        return 1;
    }
    inputFile >> numCities >> numEdges;

    // Create a graph represented by an adjacency matrix
    vector<vector<int>> graph(numCities, vector<int>(numCities, INF));

    // Initialize diagonal elements to 0 (distance from a city to itself is 0)
    for (int i = 0; i < numCities; ++i) {
        graph[i][i] = 0;
    }

    // Read and process each edge
    for (int i = 0; i < numEdges; ++i) {
        int start, end, distance;
        inputFile >> start >> end >> distance;
        graph[start][end] = distance;
        graph[end][start] = distance;
    }

    // Close the file
    inputFile.close();

    // Print the adjacency matrix
    cout << "Adjacency Matrix:" << endl;
    printGraph(graph);

    // Genetic algorithm

    return 0;
}