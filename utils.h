#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <unordered_map>
#include <limits>

using namespace std;

#define INF numeric_limits<int>::max()

struct vector_hash {
    size_t operator()(const std::vector<int>& v) const {
        size_t hash_value = 0;
        for (int element : v) {
            hash_value ^= std::hash<int>()(element) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
        }
        return hash_value;
    }
};


void printGraph(vector<vector<double>>& graph);
vector<int> kruskalMST(vector<vector<double>>& graph, vector<vector<double>>& edges, int startCity, int numCities);
vector<int> generateRandomTour(int numCities);
vector<vector<int>> initializePopulation(int N, int numCities);
void geneticAlgorithm(vector<vector<double>>& graph, vector<vector<double>>& edges, int startCity, int numCities);
vector<pair<int, int>> selection(vector<vector<int>>& population);
int selectParent(vector<vector<int>>& population);
void mutation(vector<int>& offspring, double rate);
void applyMutation(vector<int>& tour);
pair<int, int> randomTwoDifferentIndices(int size);
pair<vector<int>, vector<int>> crossover(vector<int>& parent1, vector<int>& parent2, double rate);
double evaluateFitness(vector<int>& tour, vector<vector<double>>& graph);
double calculateTourDistance(vector<int>& tour, vector<vector<double>>& graph);
void calculateFitness(vector<vector<int>>& population, vector<std::vector<double>>& graph, unordered_map<vector<int>, int, vector_hash>& fitnessMap);
vector<vector<int>> selectBestIndividuals(vector<vector<int>>& population, int count, vector<vector<double>>& graph);

#endif 
