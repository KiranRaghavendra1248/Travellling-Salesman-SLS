#ifndef UTILS_H
#define UTILS_H
#include <vector>
using namespace std;

#define INF numeric_limits<int>::max()

// Utility function declarations
void printGraph(vector<vector<double>> &graph);
vector<int> kruskalMST(vector<vector<double>> &graph, vector<vector<double>> &edges, int startCity, int numCities);
vector<int> generateRandomTour(int numCities);
vector<vector<int>> initializePopulation(int N, int numCities);
void geneticAlgorithm(vector<vector<double>>& graph,  vector<vector<double>>& edges, int startCity, int numCities);
vector<pair<int, int>> selection(vector<vector<int>>& population);
int selectParent(vector<vector<int>>& population);
void Mutation(vector<int>& offspring, double rate);
void ApplyMutation(vector<int>& tour);
pair<int, int> RandomTwoDifferentIndices(int size);
pair<vector<int>, vector<int>> Crossover(vector<int>& parent1, vector<int>& parent2, double rate);
double EvaluateFitness(vector<int>& tour, vector<vector<double>>& graph);
double CalculateTourDistance(vector<int>& tour, vector<vector<double>>& graph);
vector<vector<int>> select_best_individuals(vector<vector<int>>& population, int count, vector<vector<double>>& graph);

#endif // UTILS_H
