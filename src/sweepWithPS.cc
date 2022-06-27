// sequetial code with sweep
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <queue>
#include <climits>

using namespace std;

namespace globalTimer {
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::seconds timeLimit;
    std::chrono::high_resolution_clock::time_point endTime;

    void setLimit(std::chrono::seconds limit) {
        timeLimit = limit;
    }

    void startTimer(std::chrono::seconds limit) {
        startTime = std::chrono::high_resolution_clock::now();
        setLimit(limit);
    }

    template <class Duration = std::chrono::nanoseconds>
    Duration getTime()
    {
        endTime = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<Duration>(endTime - startTime);
    }

    bool TLE() {
        return getTime() >= timeLimit;
    }
};


struct Grid {
    int cost;
    int prev;
    Grid() : cost(INT_MAX), prev(-1) {}
};

struct Pin {
    int x;
    int y;
    bool routed;
    Pin(int _x, int _y, int r) : x(_x), y(_y), routed(r) {}
};

int W, H, numOfPins;
vector<vector<int>> hWeight, vWeight;
vector<vector<Grid *>> gridMap;
vector<Pin *> pins;
bool changed;

void parser(string filename) {
    ifstream fin(filename);
    string str;
    
    // W, H
    fin >> W >> H;

    // # of pins
    fin >> numOfPins;

    // pins
    for (int i = 0; i < numOfPins; i++) {
        fin >> str;
        if (str != "Pin") cout << "Input: Pin error" << endl;

        int x, y;
        fin >> x >> y;
        Pin *p = new Pin(x, y, false);
        pins.push_back(p);
    }

    // vertical line weight
    vWeight.resize(H);
    for (int i = 0; i < H; i++) vWeight[i].resize(W-1);

    for (int i = 0; i < H; i++) {
        fin >> str;
        if (str != "Vertical") cout << "Input: Vertical error" << endl;

        for (int j = 0; j < W - 1; j++) {
            fin >> vWeight[i][j];
        }
    }

    // horizontal line weight
    hWeight.resize(H-1);
    for (int i = 0; i < H - 1; i++) hWeight[i].resize(W);

    for (int i = 0; i < H - 1; i++) {
        fin >> str;
        if (str != "Horizontal") cout << "Input: Horizontal error" << endl;

        for (int j = 0; j < W; j++) {
            fin >> hWeight[i][j];
        }
    }
    
    // grid map
    gridMap.resize(H);
    for (int i = 0; i < H; i++) {
        gridMap[i].resize(W);
    }
    
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            gridMap[i][j] = new Grid();
        }
    }

    // for testing
    /*
    cout << "Pins\n";
    for (int i = 0; i < numOfPins; i++) {
        cout << pins[i]->x << " " << pins[i]->y << endl;
    }
    cout << endl;

    cout << "vWeight\n";
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W - 1; j++) {
            cout << vWeight[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "hWeight\n";
    for (int i = 0; i < H - 1; i++) {
        for (int j = 0; j < W; j++) {
            cout << hWeight[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    */
}

void initialize() {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i][j]->cost != 0) {
                gridMap[i][j]->cost = INT_MAX;
            }
        }
    }
}

void writeOutput(string filename) {
    ofstream fout(filename);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i][j]->cost == 0) {
                fout << i << " " << j << "\n";
            }
        }
    }
    fout.close();
}

int transTo1D(int x, int y) {
    return x * W + y;
}

void horizontalSweep() {
    for (int x = 0; x < H; x++) {
        vector<int> prefixSum(W, 0);
        for (int i = 1; i < W; i++) prefixSum[i] = prefixSum[i-1] + vWeight[x][i - 1];

        vector<int> tmp(W);
        for (int i = 0; i < W; i++) tmp[i] = gridMap[x][i]->cost - prefixSum[i];

        for (int i = 1; i < W; i++) {
            if (tmp[i] > tmp[i - 1]) {
                tmp[i] = tmp[i - 1];
                gridMap[x][i]->prev = transTo1D(x, i - 1);
                changed = true;
            }
        }

        for (int i = 0; i < W; i++) {
            gridMap[x][i]->cost = tmp[i] + prefixSum[i];
        }
    }
    for (int x = 0; x < H; x++) {
        vector<int> prefixSum(W, 0);
        for (int i = W - 2; i >= 0; i--) prefixSum[i] = prefixSum[i+1] + vWeight[x][i];

        vector<int> tmp(W);
        for (int i = 0; i < W; i++) tmp[i] = gridMap[x][i]->cost - prefixSum[i];

        for (int i = W - 2; i >= 0; i--) {
            if (tmp[i] > tmp[i + 1]) {
                tmp[i] = tmp[i + 1];
                gridMap[x][i]->prev = transTo1D(x, i + 1);
                changed = true;
            }
        }

        for (int i = 0; i < W; i++) {
            gridMap[x][i]->cost = tmp[i] + prefixSum[i];
        }
    }
}

void verticalSweep() {
    for (int y = 0; y < W; y++) {
        vector<int> prefixSum(H, 0);
        for (int i = 1; i < H; i++) prefixSum[i] = prefixSum[i-1] + hWeight[i - 1][y];

        vector<int> tmp(H);
        for (int i = 0; i < H; i++) tmp[i] = gridMap[i][y]->cost - prefixSum[i];


        for (int i = 1; i < H; i++) {
            if (tmp[i] > tmp[i - 1]) {
                tmp[i] = tmp[i - 1];
                gridMap[i][y]->prev = transTo1D(i - 1, y);
                changed = true;
            }
        }
        
        for (int i = 0; i < H; i++) {
            gridMap[i][y]->cost = tmp[i] + prefixSum[i];
        }
    }

    for (int y = 0; y < W; y++) {
        vector<int> prefixSum(H, 0);
        for (int i = H - 2; i >= 0; i--) prefixSum[i] = prefixSum[i+1] + hWeight[i][y];

        vector<int> tmp(H);
        for (int i = 0; i < H; i++) tmp[i] = gridMap[i][y]->cost - prefixSum[i];


        for (int i = H - 2; i >= 0; i--) {
            if (tmp[i] > tmp[i + 1]) {
                tmp[i] = tmp[i + 1];
                gridMap[i][y]->prev = transTo1D(i + 1, y);
                changed = true;
            }
        }
        
        for (int i = 0; i < H; i++) {
            gridMap[i][y]->cost = tmp[i] + prefixSum[i];
        }
    }
}

void getMinCostPin() {
    int minPin = -1;
    int minCost = INT_MAX;

    for (int i = 0; i < numOfPins; i++) {
        if (!pins[i]->routed && gridMap[pins[i]->x][pins[i]->y]->cost < minCost) {
            minPin = i;
            minCost = gridMap[pins[i]->x][pins[i]->y]->cost;
        }
    }

    // retrace the path & set the cost = 0
    int x = pins[minPin]->x;
    int y = pins[minPin]->y;

    while (gridMap[x][y]->cost != 0) {
        // for debug
        if (gridMap[x][y]->prev == -1) {
            cout << "getMinCost: cannot retrace the path" << endl;
        }

        gridMap[x][y]->cost = 0;
        int p = gridMap[x][y]->prev;
        x = p / W;
        y = p % W;
    }
    pins[minPin]->routed = true;

    // for testing 
    /*
    cout << "Pin: " << minPin << endl;
    cout << "Coord: " << pins[minPin]->x << " " << pins[minPin]->y << endl;
    cout << endl;
    */
}

// for testing
void printGridMap() {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i][j]->cost == INT_MAX) 
                cout << setw(3) << "INF" << " ";
            else 
                cout << setw(3) << gridMap[i][j]->cost << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void sweep() {
    int iteration = 0;
    while (true) {
        iteration++;
        changed = false;
        horizontalSweep();
        verticalSweep();
        
        if (!changed) break;
    }

    cout << "Iteration: " << iteration << endl;
}

int main(int argc, char **argv) {
    // check the argument is enough
    if (argc != 3) {
        cout << "The number of arguments are wrong!!" << endl;
        return 0;
    }

    // parse the input file
    parser(argv[1]);

    // routing (sweep)
    gridMap[pins[0]->x][pins[0]->y]->cost = 0;
    pins[0]->routed = true;
    // for testing 
    printGridMap();

    for (int i = 0; i < numOfPins - 1; i++) {
        initialize();
        sweep();
        getMinCostPin();

        // for testing 
        printGridMap();
    }

    writeOutput(argv[2]);

    return 0;
}