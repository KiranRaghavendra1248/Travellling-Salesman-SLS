#include <iostream>
#include <limits>
#include <set>
#include <vector>
#include <iomanip>
#include <cassert>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include "utils.h"

using namespace std;

// Utility classes and struct

struct VectorComparator{
    bool operator()(const vector<int> &a, const vector<int> &b) const
    {
        return a < b; // Use lexicographical comparison for vectors
    }
};

class DisjointSet{
private:
    vector<int> parent, rank;

public:
    DisjointSet(int n){
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i){
            parent[i] = i;
        }
    }

    // Find the set (representative) of the given element
    int findSet(int x){
        if (x != parent[x]){
            parent[x] = findSet(parent[x]);
        }
        return parent[x];
    }

    // Merge two sets if they are different
    void unionSets(int x, int y){
        int rootX = findSet(x);
        int rootY = findSet(y);

        if (rootX != rootY)
        {
            if (rank[rootX] < rank[rootY]){
                swap(rootX, rootY);
            }
            parent[rootY] = rootX;
            if (rank[rootX] == rank[rootY]){
                rank[rootX]++;
            }
        }
    }
};

// Utility function definitions

void printGraph(vector<vector<double>> &graph){
    for (auto &row : graph){
        for (double distance : row){
            if (distance == INF){
                cout << setw(5) << "INF";
            }
            else{
                cout << setw(5) << distance;
            }
        }
        cout << endl;
    }
}

vector<int> kruskalMST(vector<vector<double>> &graph, vector<vector<double>> &edges, int startCity, int numCities){
    int n = numCities;

    // Create disjoint set for each city
    DisjointSet ds(n);

    // Kruskal's algorithm to find MST
    vector<vector<double>> mstEdges;
    for (auto &edge : edges){
        int u = edge[0];
        int v = edge[1];

        // Check if adding this edge forms a cycle in the MST
        if (ds.findSet(u) != ds.findSet(v)){
            // Add the edge to the MST
            mstEdges.push_back(edge);
            // Merge the sets of vertices u and v
            ds.unionSets(u, v);
        }
    }

    // Create an adjacency list to represent the MST
    vector<vector<double>> adjacencyList(n);
    for (auto &edge : mstEdges){
        adjacencyList[edge[0]].push_back(edge[1]);
        adjacencyList[edge[1]].push_back(edge[0]);
    }

    // Perform a depth-first traversal to create a graph path
    vector<int> graphPath;
    vector<bool> visited(n, false);

    // Depth-first traversal function
    function<void(int)> dfs = [&](int current){
        graphPath.push_back(current);
        visited[current] = true;
        for (double neighbor : adjacencyList[current]){
            if (!visited[neighbor]){
                dfs(neighbor);
            }
        }
    };
    dfs(startCity);
    return graphPath;
}

vector<int> generateRandomTour(int numCities){
    vector<int> tour(numCities);
    for (int i = 0; i < numCities; ++i){
        tour[i] = i;
    }
    // Use a random device as a seed for the random number generator
    random_device rd;
    // Use the random device to seed the random number generator
    mt19937 gen(rd());
    // Shuffle the vector to generate a random permutation
    shuffle(tour.begin()+1, tour.end(), gen);
    return tour;
}

vector<vector<int>> initializePopulation(int N, int numCities){
    cout << "Initializing population" << endl;
    vector<vector<int>> generatedPopulation;
    for (int i = 0; i < N; ++i){
        vector<int> vec = generateRandomTour(numCities);
        generatedPopulation.push_back(vec);
    }
    // Generate vector from set
    return generatedPopulation;
}

vector<pair<int,int>>selection(vector<vector<int>> &population){
    vector<pair<int, int>> selectedParents;
    for (size_t i = 0; i < population.size(); ++i){
        int parent1 = selectParent(population);
        int parent2 = selectParent(population);
        selectedParents.push_back(make_pair(parent1, parent2));
    }
    return selectedParents;
}

int selectParent(vector<vector<int>> &population){
    int randomIndex = rand() % population.size();
    return randomIndex;
}

void Mutation(vector<int> &offspring, double rate){
    int randMax = RAND_MAX + 1;
    int tourLength = offspring.size(); // Assuming offspring is a vector of integers representing tours
    for (int i = 0; i < tourLength; ++i){
        int randomNumber = rand();
        double probability = static_cast<double>(randomNumber) / randMax;
        if (probability < rate)
        {
            ApplyMutation(offspring);
        }
    }
}

void ApplyMutation(vector<int> &tour){
    int tourLength = tour.size();
    pair<int, int> indices = RandomTwoDifferentIndices(tourLength);
    swap(tour[indices.first], tour[indices.second]);
}

pair<int, int> RandomTwoDifferentIndices(int size){
    int index1 = (rand() % (size-1)) + 1;
    int index2;
    do{
        index2 = (rand() % (size-1)) + 1;
    } while (index2 == index1);
    return make_pair(index1, index2);
}

pair<vector<int>,vector<int>> Crossover(vector<int>& parent1, vector<int>& parent2, double rate){
    if ((rand() % 100) < rate * 100) {
        assert (parent1.size() == parent2.size());
        size_t crossoverPoint = (rand() % (min(parent1.size(), parent2.size())-1)) + 1;
        vector<int> offspring1(parent1.size());
        vector<int> offspring2(parent2.size());
        unordered_set<int> elementsInOffspring1;
        for (size_t i = 0; i < crossoverPoint; ++i) {
            offspring1[i] = parent1[i];
            elementsInOffspring1.insert(parent1[i]);
        }
        size_t offspring1Index = crossoverPoint;
        for (size_t i = 0; i < parent2.size(); ++i) {
            int element = parent2[i];
            if (elementsInOffspring1.find(element) == elementsInOffspring1.end()) {
                offspring1[offspring1Index++] = element;
                elementsInOffspring1.insert(element);
            }
        }
        unordered_set<int> elementsInOffspring2;
        for (size_t i = 0; i < crossoverPoint; ++i) {
            offspring2[i] = parent2[i];
            elementsInOffspring2.insert(parent2[i]);
        }
        size_t offspring2Index = crossoverPoint;
        for (size_t i = 0; i < parent1.size(); ++i) {
            int element = parent1[i];
            if (elementsInOffspring2.find(element) == elementsInOffspring2.end()) {
                offspring2[offspring2Index++] = element;
                elementsInOffspring2.insert(element);
            }
        }
        return make_pair(offspring1, offspring2);
    } 
    else {
        return make_pair(parent1, parent2);
    }
}

double CalculateTourDistance(vector<int> &tour, vector<vector<double>> &graph){
    double total_distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
    {
        int current_city = tour[i];
        int next_city = tour[i + 1];
        double distance = graph[current_city][next_city];
        total_distance += distance;
    }
    int last_city = tour.back();
    int first_city = tour.front();
    double distance_to_start = graph[last_city][first_city];
    total_distance += distance_to_start;
    return total_distance;
}

double EvaluateFitness(vector<int> &tour, vector<vector<double>> &graph){
    double distance = CalculateTourDistance(tour, graph);
    double fitness = (distance > 0) ? (1.0 / distance) : numeric_limits<double>::infinity();
    return fitness;
}

vector<vector<int>> select_best_individuals(vector<vector<int>> &population, int count, vector<vector<double>> &graph){
    sort(population.begin(), population.end(), [&](vector<int> &a, vector<int> &b)
         {return EvaluateFitness(a, graph) < EvaluateFitness(b, graph); });
    if (count <= population.size()){
        return vector<vector<int>>(population.begin(), population.begin() + count);
    }
    else{
        return population;
    }
}

void CalculateFitness(vector<vector<int>> &population, vector<vector<double>> &graph, unordered_map<vector<int>, int, vector_hash>& fitnessMap){
    for (vector<int>& tour : population) {
        int fitness = EvaluateFitness(tour, graph);
        // Store the fitness in the hashmap with the tour as the key
        fitnessMap[tour] = fitness;
    }
}

vector<vector<int>> SelectBestSolution(vector<vector<int>> &population, int count, vector<vector<double>> &graph){
    // Declare an unordered_map to store tour-fitness pairs
    unordered_map<vector<int>, int, vector_hash> fitnessMap;

    // Calculate fitness for each individual in the population
    CalculateFitness(population, graph, fitnessMap);

    // Sort the population based on fitness (lower fitness values come first)
    sort(population.begin(), population.end(), [&fitnessMap](vector<int>& a,vector<int>& b){
        return fitnessMap[a] < fitnessMap[b];
    });

    if (count <= population.size()){
        return vector<vector<int>>(population.begin(), population.begin() + count);
    }
    else{
        return population; // Return the whole population if count exceeds population size
    }
}

void geneticAlgorithm(vector<vector<double>> &graph, vector<vector<double>> &edges, int startCity, int numCities){
    int populationSize;
    vector<int> initialSolution;
    vector<vector<int>> population, generatedRandomizedPopulation;

    // Generate inital solution using Krustal's MST algo
    initialSolution = kruskalMST(graph, edges, startCity, numCities);

    // Initalize population
    population.push_back(initialSolution);
    generatedRandomizedPopulation = initializePopulation(populationSize - 1, numCities);
    population.insert(population.end(), generatedRandomizedPopulation.begin(), generatedRandomizedPopulation.end());

    // Params
    int count = 1;
    populationSize = 25;
    int generations = 25;
    double crossoverRate = 0.4;
    double mutationRate = 0.2;

    for (int generation = 0; generation < generations; ++generation){
        cout << "Running generation " << generation << endl;
        vector<pair<int, int>> selectedParents = selection(population);
        vector<vector<int>> offspring;
        for (size_t i = 0; i < selectedParents.size(); ++i){
            pair<int, int> parents = selectedParents[i];
            vector<int> &parent1 = population[parents.first];
            vector<int> &parent2 = population[parents.second];
            pair<vector<int>, vector<int>> child = Crossover(parent1, parent2, crossoverRate);
            Mutation(child.first, mutationRate);
            Mutation(child.second, mutationRate);
            offspring.push_back(child.first);
            offspring.push_back(child.second);
        }

        population.insert(population.end(), offspring.begin(), offspring.end());
        population = select_best_individuals(population, populationSize, graph);
    }

    vector<vector<int>> finalSolution = SelectBestSolution(population, count, graph);
    cout << "Final Travelling Salesman Solution:" << endl;
    for (const auto &tour : finalSolution){
        for (int city : tour){
            cout << city << " ";
        }
        cout << startCity << " ";
        cout << endl;
    }
}
