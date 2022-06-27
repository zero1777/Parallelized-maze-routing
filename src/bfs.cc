// sequetial code with bfs
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

    void startTimer() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    template <class Duration = std::chrono::nanoseconds>
    Duration getTime()
    {
        endTime = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<Duration>(endTime - startTime);
    }

    template <class Duration = std::chrono::nanoseconds>
    Duration endTimer() {
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

void bfs() {
    queue<int> qu;    

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i][j]->cost == 0) {
                qu.push(transTo1D(i, j));
            }
        }
    }

    while (!qu.empty()) {
        int q = qu.front();
        qu.pop();

        int x = q / W;
        int y = q % W;

        // right 
        if (y < W - 1) {
            if (gridMap[x][y + 1]->cost > gridMap[x][y]->cost + vWeight[x][y]) {
                gridMap[x][y + 1]->cost = gridMap[x][y]->cost + vWeight[x][y];
                gridMap[x][y + 1]->prev = transTo1D(x, y);
                qu.push(transTo1D(x, y + 1));
            }
        }

        // left 
        if (y > 0) {
            if (gridMap[x][y - 1]->cost > gridMap[x][y]->cost + vWeight[x][y - 1]) {
                gridMap[x][y - 1]->cost = gridMap[x][y]->cost + vWeight[x][y - 1];
                gridMap[x][y - 1]->prev = transTo1D(x, y);
                qu.push(transTo1D(x, y - 1));
            }
        }

        // down
        if (x < H - 1) {
            if (gridMap[x + 1][y]->cost > gridMap[x][y]->cost + hWeight[x][y]) {
                gridMap[x + 1][y]->cost = gridMap[x][y]->cost + hWeight[x][y];
                gridMap[x + 1][y]->prev = transTo1D(x, y);
                qu.push(transTo1D(x + 1, y));
            }
        }

        // up 
        if (x > 0) {
            if (gridMap[x - 1][y]->cost > gridMap[x][y]->cost + hWeight[x - 1][y]) {
                gridMap[x - 1][y]->cost = gridMap[x][y]->cost + hWeight[x - 1][y];
                gridMap[x - 1][y]->prev = transTo1D(x, y);
                qu.push(transTo1D(x - 1, y));
            }
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

int main(int argc, char **argv) {
    // check the argument is enough
    if (argc != 3) {
        cout << "The number of arguments are wrong!!" << endl;
        return 0;
    }

    globalTimer::startTimer();
    // parse the input file
    parser(argv[1]);
    std::chrono::nanoseconds _input = globalTimer::endTimer();

    // routing (bfs)
    gridMap[pins[0]->x][pins[0]->y]->cost = 0;
    pins[0]->routed = true;
    // for testing 
    // printGridMap();

    globalTimer::startTimer();
    for (int i = 0; i < numOfPins - 1; i++) {
        initialize();
        bfs();
        getMinCostPin();

        // for testing 
        // printGridMap();
    }
    std::chrono::nanoseconds _kernel = globalTimer::endTimer();
    printGridMap();

    globalTimer::startTimer();
    writeOutput(argv[2]);
    std::chrono::nanoseconds _output = globalTimer::endTimer();

    // cout time usage
    auto inputTime = _input.count();
    auto outputTime = _output.count();
    auto kernelTime = _kernel.count();

    cout << "I/O time: " << inputTime / 1e9 + outputTime / 1e9 << endl;
    cout << "Kernel time: " << kernelTime / 1e9 << endl;

    return 0;
}