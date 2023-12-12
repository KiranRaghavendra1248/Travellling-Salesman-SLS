#include <iostream> 
#include <limits> 
#include <vector>
#include <iomanip> 
#include "utils.h" 

using namespace std;

int INF = numeric_limits<int>::max();

// Utility functions definition

void printGraph(const vector<vector<int>>& graph){
    for (const auto& row : graph) {
        for (int distance : row) {
            if (distance == INF) {
                cout << setw(5) << "INF";
            } else {
                cout << setw(5) << distance;
            }
        }
        cout << endl;
    }
}