#ifndef UTILS_H
#define UTILS_H
#include <vector>
using namespace std;

#define INF numeric_limits<int>::max()

// Utility function declarations
void printGraph(vector<vector<int>>& graph);
vector<int> kruskalMST(vector<vector<int>>& graph,  vector<vector<int>>& edges, int startCity, int numCities, int numEdges);
vector<int> generateRandomTour(int numCities);
vector<vector<int>> initializePopulation(int N, int numCities);
void geneticAlgorithm(vector<vector<int>>& graph,  vector<vector<int>>& edges, int startCity, int numCities, int numEdges);

#endif // UTILS_H