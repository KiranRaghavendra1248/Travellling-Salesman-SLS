#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cassert>
#include <iomanip>
#include "utils.h"

using namespace std;

#define INF numeric_limits<int>::max()

int main(int argc, char* argv[]) {
    int numCities, startCity;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    const char* filename = argv[1];
    ifstream inputFile(filename);

    if (!inputFile) {
        cerr << "Error: Unable to open file '" << filename << "'" << std::endl;
        return 1;
    }
    inputFile >> numCities;

    vector<vector<double>> graph(numCities, vector<double>(numCities, INF));
    vector<vector<double>> edges;

    for (double i = 0; i < numCities; i++) {
        for (double j = 0; j < numCities; j++) {
            inputFile >> graph[i][j];
            edges.push_back({i, j, graph[i][j]});
            edges.push_back({j, i, graph[i][j]});
        }
    }

    inputFile.close();

    cout << "Adjacency Matrix:" << endl;
    printGraph(graph);

    startCity = 0;
    clock_t start_time = clock();
    geneticAlgorithm(graph, edges, startCity, numCities);
    clock_t end_time = clock();

    double t = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Time Taken : "<< t << " seconds" << endl;

    
    return 0;
}
