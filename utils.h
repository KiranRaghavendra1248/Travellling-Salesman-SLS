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
void geneticAlgorithm(vector<vector<int>>& graph,  vector<vector<int>>& edges, int startCity, int numCities, int numEdges, vector<vector<int>>& TSP);
vector<pair<int, int>> selection(vector<vector<int>>& population);
int selectParent(vector<vector<int>>& population);
void Mutation(vector<int>& offspring, double rate);
void ApplyMutation(vector<int>& tour);
pair<int, int> RandomTwoDifferentIndices(int size);
pair<vector<int>, vector<int>> Crossover(vector<int>& parent1, vector<int>& parent2, double rate);
double EvaluateFitness(vector<int>& tour, vector<std::vector<int>>& TSP);
int CalculateTourDistance(vector<int>& tour, vector<std::vector<int>>& TSP);
vector<vector<int>> select_best_individuals(vector<std::vector<int>>& population, int count, vector<std::vector<int>>& TSP);

#endif // UTILS_H
