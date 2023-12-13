#include <iostream>
#include <fstream> 
#include <vector>
#include <limits> 
#include <cassert>
#include <iomanip> 
#include "utils.h"
using namespace std;

#define INF numeric_limits<int>::max()
 
int main(int argc, char* argv[])
{   
    int numCities, startCity;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1; // Exit with an error code
    }

    // Open the file specified in the command-line argument
    const char* filename = argv[1];
    ifstream inputFile(filename);

    // Check if the file is successfully opened
    if (!inputFile) {
        cerr << "Error: Unable to open file '" << filename << "'" << std::endl;
        return 1; // Exit with an error code
    }
    inputFile >> numCities;


    // Create a graph represented by an adjacency matrix
    vector<vector<double>> graph(numCities, vector<double>(numCities, INF));
    vector<vector<double>> edges;

    for (double i = 0; i < numCities; i++) {
        for (double j = 0; j < numCities; j++) {
            inputFile >> graph[i][j];
            edges.push_back({i, j, graph[i][j]});
            edges.push_back({j, i, graph[i][j]});
        }
    }

    // Close the file
    inputFile.close();

    // Print the adjacency matrix
    cout << "Adjacency Matrix:" << endl;
    printGraph(graph);

    // Let start_city = 0 by default
    startCity = 0;

    // Genetic algorithm
    geneticAlgorithm(graph, edges, startCity, numCities);

    return 0;
}
