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

struct VectorComparator {
    bool operator()(const vector<int> &a, const vector<int> &b) const {
        return a < b;
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

    int findSet(int x) {
        if (x != parent[x]) {
            parent[x] = findSet(parent[x]);
        }
        return parent[x];
    }

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

void printGraph(vector<vector<double>> &graph) {
    for (auto &row : graph) {
        for (double distance : row) {
            if (distance == INF) {
                cout << setw(5) << "INF";
            } else {
                cout << setw(5) << distance;
            }
        }
        cout << endl;
    }
}

vector<int> kruskalMST(vector<vector<double>> &graph, vector<vector<double>> &edges, int startCity, int numCities) {
    int n = numCities;
    DisjointSet ds(n);
    vector<vector<double>> mstEdges;

    for (auto &edge : edges) {
        int u = edge[0];
        int v = edge[1];

        if (ds.findSet(u) != ds.findSet(v)) {
            mstEdges.push_back(edge);
            ds.unionSets(u, v);
        }
    }

    vector<vector<double>> adjacencyList(n);
    for (auto &edge : mstEdges) {
        adjacencyList[edge[0]].push_back(edge[1]);
        adjacencyList[edge[1]].push_back(edge[0]);
    }

    vector<int> graphPath;
    vector<bool> visited(n, false);

    function<void(int)> dfs = [&](int current) {
        graphPath.push_back(current);
        visited[current] = true;
        for (double neighbor : adjacencyList[current]) {
            if (!visited[neighbor]) {
                dfs(neighbor);
            }
        }
    };
    dfs(startCity);
    return graphPath;
}

vector<int> generateRandomTour(int numCities) {
    vector<int> tour(numCities);
    for (int i = 0; i < numCities; ++i) {
        tour[i] = i;
    }
    random_device rd;
    mt19937 gen(rd());
    shuffle(tour.begin() + 1, tour.end(), gen);
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

/*vector<vector<int>> initializePopulation(int N, int numCities) {
    cout << "Initializing population" << endl;
    vector<vector<int>> generatedPopulation;
    for (int i = 0; i < N; ++i) {
        vector<int> vec = generateRandomTour(numCities);
        generatedPopulation.push_back(vec);
    }
    return generatedPopulation;
}*/

vector<pair<int, int>> selection(vector<vector<int>> &population) {
    vector<pair<int, int>> selectedParents;
    for (size_t i = 0; i < population.size(); ++i) {
        int parent1 = selectParent(population);
        int parent2;
        do {
            parent2 = selectParent(population);
        } while (parent2 == parent1);
        selectedParents.push_back(make_pair(parent1, parent2));
    }
    return selectedParents;
}

int selectParent(vector<vector<int>> &population) {
    int randomIndex = rand() % population.size();
    return randomIndex;
}

void mutation(vector<int> &offspring, double rate) {
    int randMax = RAND_MAX + 1;
    int tourLength = offspring.size();
    for (int i = 0; i < tourLength; ++i) {
        int randomNumber = rand();
        double probability = static_cast<double>(randomNumber) / randMax;
        if (probability < rate) {
            applyMutation(offspring);
           /* cout << "\nAfter Mutation: ";
            for (int city : offspring) {
                cout << city << " ";
            }
           cout << endl; */
        }
    }
}

void applyMutation(vector<int> &tour) {
    int tourLength = tour.size();
    pair<int, int> indices = randomTwoDifferentIndices(tourLength);
    swap(tour[indices.first], tour[indices.second]);
}

pair<int, int> randomTwoDifferentIndices(int size) {
    int index1 = (rand() % (size - 1)) + 1;
    int index2 = (rand() % (size - 1)) + 1;
    return make_pair(index1, index2);
}

pair<vector<int>, vector<int>> crossover(vector<int> &parent1, vector<int> &parent2, double rate) {
    if ((rand() % 100) < rate * 100) {
        assert(parent1.size() == parent2.size());
        size_t crossoverPoint = (rand() % (min(parent1.size(), parent2.size()) - 1)) + 1;
       /* cout << "Parent 1: ";
        for (int city : parent1) {
            cout << city << " ";
        }
        cout << "\nParent 2: ";
        for (int city : parent2) {
            cout << city << " ";
        }*/
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
       /* cout << "\nOffspring 1: ";
        for (int city : offspring1) {
            cout << city << " ";
        }
        cout << "\nOffspring 2: ";
        for (int city : offspring2) {
            cout << city << " ";
        }
        cout << endl;*/
        return make_pair(offspring1, offspring2);
    } else {
        return make_pair(parent1, parent2);
    }
}



double calculateTourDistance(vector<int> &tour, vector<vector<double>> &graph) {
    double total_distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
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

double evaluateFitness(vector<int> &tour, vector<vector<double>> &graph) {
    double distance = calculateTourDistance(tour, graph);
    // cout << "Minimum distance of the tour : " << calculateTourDistance(tour, graph) << endl;
    double fitness = (distance > 0) ? (1.0 / distance) : numeric_limits<double>::infinity();
    return fitness;
}

vector<vector<int>> selectBestIndividuals(vector<vector<int>> &population, int count, vector<vector<double>> &graph) {
    sort(population.begin(), population.end(), [&](vector<int> &a, vector<int> &b) {
        //return evaluateFitness(a, graph) < evaluateFitness(b, graph);
        return calculateTourDistance(a, graph) < calculateTourDistance(b, graph);
    });
    if (count <= population.size()) {
        return vector<vector<int>>(population.begin(), population.begin() + count);
    } else {
        return population;
    }
}

vector<int> selectBestSolution(vector<vector<int>> &population, vector<vector<double>> &graph) {
    sort(population.begin(), population.end(), [&](vector<int> &a, vector<int> &b) {
        return calculateTourDistance(a, graph) < calculateTourDistance(b, graph);
    });
    return population.front();
}

void geneticAlgorithm(vector<vector<double>> &graph, vector<vector<double>> &edges, int startCity, int numCities) {
    vector<int> initialSolution;
    vector<vector<int>> population, generatedRandomizedPopulation;

    int populationSize = 250;
    
    int generations = 1000;
    double crossoverRate = 1;
    double mutationRate = 0.6;

    initialSolution = kruskalMST(graph, edges, startCity, numCities);

    population.push_back(initialSolution);
    generatedRandomizedPopulation = initializePopulation(populationSize - 1, numCities);
    population.insert(population.end(), generatedRandomizedPopulation.begin(), generatedRandomizedPopulation.end());

    for (int generation = 0; generation < generations; ++generation) {
         cout << "Running generation " << generation << endl;
        vector<pair<int, int>> selectedParents = selection(population);
        vector<vector<int>> offspring;
        for (size_t i = 0; i < selectedParents.size(); ++i) {
            pair<int, int> parents = selectedParents[i];
            vector<int> &parent1 = population[parents.first];
            vector<int> &parent2 = population[parents.second];
            pair<vector<int>, vector<int>> child = crossover(parent1, parent2, crossoverRate);
            mutation(child.first, mutationRate);
            mutation(child.second, mutationRate);
            offspring.push_back(child.first);
            offspring.push_back(child.second);
        }
        population.insert(population.end(), offspring.begin(), offspring.end());
        population = selectBestIndividuals(population, populationSize, graph);
    }

    vector<int> finalSolution = selectBestSolution(population, graph);
    cout << "Final Travelling Salesman Solution:" << endl;
    for (auto &tour : finalSolution) {
            cout << tour << " ";
        }
        cout << startCity << " ";
        cout << endl;
        cout << "Minimum distance of the tour : " << calculateTourDistance(finalSolution, graph) << endl;
        
    }

