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

#define sharedSize 1024
#define NumOfThreads 128

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

void parse_basic(string filename, int &W, int &H, int &numOfPins) {
    ifstream fin(filename);

    // W, H
    fin >> W >> H;

    // # of pins
    fin >> numOfPins;
}

void parser(string filename, int W, int H, int numOfPins, int3 *pins, int *hWeight, int *vWeight, int2 *gridMap) {
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
        pins[i].x = x;
        pins[i].y = y;
        pins[i].z = 0;
    }

    // vertical line weight
    for (int i = 0; i < H; i++) {
        fin >> str;
        if (str != "Vertical") cout << "Input: Vertical error" << endl;

        for (int j = 0; j < W - 1; j++) {
            fin >> vWeight[i * W + j];
        }
        vWeight[i * W + W - 1] = 0;
    }

    // horizontal line weight
    for (int i = 0; i < H - 1; i++) {
        fin >> str;
        if (str != "Horizontal") cout << "Input: Horizontal error" << endl;

        for (int j = 0; j < W; j++) {
            fin >> hWeight[i * W + j];
        }
    }
    // H
    for (int j = 0; j < W; j++) {
        hWeight[(H - 1) * W + j] = 0;
    }

    // grid map
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            gridMap[i * W + j].x = INT_MAX;
            gridMap[i * W + j].y = -1;
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
        for (int j = 0; j < W; j++) {
            cout << vWeight[i * W + j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "hWeight\n";
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            cout << hWeight[i * W + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    */
}

__global__ void initialize(int H, int W, int2 *gridMap) {
    int tid = threadIdx.x;
    int stride = blockDim.x;
    int i = blockIdx.x;

    if (i >= H) return ;

    for (int j = tid; j < W; j += stride) {
        if (gridMap[i * W + j].x != 0) {
            gridMap[i * W + j].x = INT_MAX;
        }
    }
}

void writeOutput(string filename, int H, int W, int2 *gridMap) {
    ofstream fout(filename);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i * W + j].x == 0) {
                fout << i << " " << j << "\n";
            }
        }
    }
    fout.close();
}

__device__ __host__ int transTo1D(int x, int y, int W) {
    return x * W + y;
}

__global__ void horizontalSweep1(int H, int W, int2 *gridMap, int *vWeight, int *changed) {
    int x = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];
    __shared__ int prefixSum[sharedSize];
    __shared__ int tmp[sharedSize];

    if (id >= W || x >= H) return ;

    for (int i = id; i < W; i += stride) {
        sharedGridMap[i] = gridMap[x * W + i];
        tmpWeight[i] = vWeight[x * W  + i];
    }

    __syncthreads();

    int offset = 1; 
    int n = sharedSize;

    for (int i = id; i < n; i += stride) {
        if (i >= W) prefixSum[i] = 0;
        else prefixSum[i] = tmpWeight[i];
    }

    __syncthreads();

	for (int d = n / 2; d > 0; d /= 2)	// build sum in place up the tree 
	{ 	
	    if (id < d) { 
        	int ai = offset * (2 * id + 1) - 1;
        	int bi = offset * (2 * id + 2) - 1;
        
            if (prefixSum[ai] != INT_MAX)
        	    prefixSum[bi] += prefixSum[ai];
        }
        
        offset *= 2;
	    __syncthreads();    
    }

    if (id == 0)
    {
    	prefixSum[n - 1] = 0;
    }	// clear the last element  

    for (int d = 1; d < n; d *= 2)	// traverse down tree &build scan 
    { 	
        offset /= 2;      
        __syncthreads();      
        if (id < d) { 

        	int ai = offset * (2 * id + 1) - 1;
            int bi = offset * (2 * id + 2) - 1;
            
            int t = prefixSum[ai];
            prefixSum[ai] = prefixSum[bi];
            if (t != INT_MAX)
                prefixSum[bi] += t;
        }
    }
    __syncthreads();

    for (int i = id; i < W; i += stride){
        tmp[i] = sharedGridMap[i].x - prefixSum[i];
    }

    __syncthreads();

    // if (id == 0) {
    //     for (int i = 1; i < W; i++) {
    //         if (tmp[i] > tmp[i - 1]) {
    //             tmp[i] = tmp[i - 1];
    //             sharedGridMap[i].y = transTo1D(x, i - 1, W);
    //             *changed = 1;
    //         }
    //     }
    // }

    int startIdx = 0;
    for (int i = id; i < W; i += stride) {
        if (i >= 1) {
            for (int s = startIdx; s < i; s++) {
                if (tmp[i] > tmp[s]) {
                    tmp[i] = tmp[s];
                    sharedGridMap[i].y = transTo1D(x, i - 1, W);
                    *changed = 1;
                }
            }
        }
        if (startIdx == 0) {
            startIdx = stride - 1;
        }
        else {
            startIdx += stride;
        }
        __syncthreads();
    }

    __syncthreads();

    for (int i = id; i < W; i += stride) {
        sharedGridMap[i].x = tmp[i] + prefixSum[i];
        gridMap[x * W + i] = sharedGridMap[i];
    }
}

__global__ void horizontalSweep2(int H, int W, int2 *gridMap, int *vWeight, int *changed) {
    int x = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];
    __shared__ int prefixSum[sharedSize];
    __shared__ int rprefixSum[sharedSize];
    __shared__ int tmp[sharedSize];

    if (id >= W || x >= H) return ;

    for (int i = id; i < W; i += stride) {
        sharedGridMap[i] = gridMap[x * W + i];
        tmpWeight[i] = vWeight[x * W  + i];
    }

    __syncthreads();

    int offset = 1; 
    int n = stride;

    for (int i = id; i < n; i += stride) {
        if (i >= W) prefixSum[i] = 0;
        else rprefixSum[i] = tmpWeight[W - i - 1];
    }

    __syncthreads();

	for (int d = n / 2; d > 0; d /= 2)	// build sum in place up the tree 
	{ 	
	    if (id < d) { 
        	int ai = offset * (2 * id + 1) - 1;
        	int bi = offset * (2 * id + 2) - 1;
        
            if (rprefixSum[ai] != INT_MAX)
        	    rprefixSum[bi] += rprefixSum[ai];
        }
        
        offset *= 2;
	    __syncthreads();    
    }

    if (id == 0)
    {
    	rprefixSum[n - 1] = 0;
    }	// clear the last element  

    for (int d = 1; d < n; d *= 2)	// traverse down tree &build scan 
    { 	
        offset /= 2;      
        __syncthreads();      
        if (id < d) { 

        	int ai = offset * (2 * id + 1) - 1;
            int bi = offset * (2 * id + 2) - 1;
            
            int t = rprefixSum[ai];
            rprefixSum[ai] = rprefixSum[bi];
            if (t != INT_MAX)
                rprefixSum[bi] += t;
        }
    }
    __syncthreads();

    for (int i = id; i < W; i += stride) {
        prefixSum[i] = rprefixSum[W - i - 1];
    }

    __syncthreads();

    for (int i = id; i < W; i += stride){
        tmp[i] = sharedGridMap[i].x - prefixSum[i];
    }

    __syncthreads();

    // if (id == 0) {
    //     for (int i = 1; i < W; i++) {
    //         if (tmp[i] > tmp[i - 1]) {
    //             tmp[i] = tmp[i - 1];
    //             sharedGridMap[i].y = transTo1D(x, i - 1, W);
    //             *changed = 1;
    //         }
    //     }
    // }

    int startIdx = W - 1;
    for (int i = id; i < W; i += stride) {
        int ci = W - i - 1;
        if (ci >= 1) {
            for (int s = startIdx; s > ci; s--) {
                if (tmp[ci] > tmp[s]) {
                    tmp[ci] = tmp[s];
                    sharedGridMap[ci].y = transTo1D(x, ci + 1, W);
                    *changed = 1;
                }
            }
        }
        if (startIdx == W - 1) {
            startIdx -= (stride - 1);
        }
        else {
            startIdx -= stride;
        }
        __syncthreads();
    }

    __syncthreads();

    for (int i = id; i < W; i += stride) {
        sharedGridMap[i].x = tmp[i] + prefixSum[i];
        gridMap[x * W + i] = sharedGridMap[i];
    }
}

__global__ void horizontalSweep3(int H, int W, int2 *gridMap, int *vWeight, int *changed) {
    int x = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];

    if (id >= W || x >= H) return ;

    for (int i = id; i < W; i += stride) {
        sharedGridMap[i] = gridMap[x * W + i];
        tmpWeight[i] = vWeight[x * W  + i];
    }

    __syncthreads();

    for (int i = 1; i <= ceil(log2f(W)); i++) {
        for (int j = id; j < W; j += stride) {
            int invert_j = W - 1 - j;
            if (invert_j % int(pow(2, i)) >= pow(2, i - 1)) {
                int idx = W - 1 - (invert_j - (invert_j % int(pow(2, i - 1))) - 1);
                int partialCost = 0;

                for (int k = idx - 1; k >= j; k--) {
                    partialCost += tmpWeight[k];
                }

                if (sharedGridMap[idx].x != INT_MAX && sharedGridMap[j].x > sharedGridMap[idx].x + partialCost) {
                    sharedGridMap[j].x = sharedGridMap[idx].x + partialCost;
                    sharedGridMap[j].y = transTo1D(x, j + 1, W);
                    *changed = 1;
                }
            }
        }
    }

    __syncthreads();

    for (int i = id; i < W; i += stride) {
        gridMap[x * W + i] = sharedGridMap[i];
    }
}

__global__ void verticalSweep1(int H, int W, int2 *gridMap, int *hWeight, int *changed) {
    int y = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];
    __shared__ int prefixSum[sharedSize];
    __shared__ int tmp[sharedSize];

    if (id >= H || y >= W) return ;

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i] = gridMap[i * W + y];
        tmpWeight[i] = hWeight[i * W + y];
    }

    __syncthreads();

    int offset = 1; 
    int n = stride;

    for (int i = id; i < n; i += stride) {
        if (i >= H) prefixSum[i] = 0;
        else prefixSum[i] = tmpWeight[i];
    }

    __syncthreads();

	for (int d = n / 2; d > 0; d /= 2)	// build sum in place up the tree 
	{ 	
	    if (id < d) { 
        	int ai = offset * (2 * id + 1) - 1;
        	int bi = offset * (2 * id + 2) - 1;
        
            if (prefixSum[ai] != INT_MAX)
        	    prefixSum[bi] += prefixSum[ai];
        }
        
        offset *= 2;
	    __syncthreads();    
    }

    if (id == 0)
    {
    	prefixSum[n - 1] = 0;
    }	// clear the last element  

    for (int d = 1; d < n; d *= 2)	// traverse down tree &build scan 
    { 	
        offset /= 2;      
        __syncthreads();      
        if (id < d) { 

        	int ai = offset * (2 * id + 1) - 1;
            int bi = offset * (2 * id + 2) - 1;
            
            int t = prefixSum[ai];
            prefixSum[ai] = prefixSum[bi];
            if (t != INT_MAX)
                prefixSum[bi] += t;
        }
    }
    __syncthreads();

    for (int i = id; i < H; i += stride){
        tmp[i] = sharedGridMap[i].x - prefixSum[i];
    }

    __syncthreads();

    // if (id == 0) {
    //     for (int i = 1; i < W; i++) {
    //         if (tmp[i] > tmp[i - 1]) {
    //             tmp[i] = tmp[i - 1];
    //             sharedGridMap[i].y = transTo1D(x, i - 1, W);
    //             *changed = 1;
    //         }
    //     }
    // }

    int startIdx = 0;
    for (int i = id; i < H; i += stride) {
        if (i >= 1) {
            for (int s = startIdx; s < i; s++) {
                if (tmp[i] > tmp[s]) {
                    tmp[i] = tmp[s];
                    sharedGridMap[i].y = transTo1D(i - 1, y, W);
                    *changed = 1;
                }
            }
        }
        if (startIdx == 0) {
            startIdx = stride - 1;
        }
        else {
            startIdx += stride;
        }
        __syncthreads();
    }

    __syncthreads();

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i].x = tmp[i] + prefixSum[i];
        gridMap[i * W + y] = sharedGridMap[i];
    }
}

__global__ void verticalSweep2(int H, int W, int2 *gridMap, int *hWeight, int *changed) {
    int y = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];
    __shared__ int prefixSum[sharedSize];
    __shared__ int rprefixSum[sharedSize];
    __shared__ int tmp[sharedSize];

    if (id >= H || y >= W) return ;

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i] = gridMap[i * W + y];
        tmpWeight[i] = hWeight[i * W  + y];
    }

    __syncthreads();

    int offset = 1; 
    int n = sharedSize;

    for (int i = id; i < n; i += stride) {
        if (i >= H) prefixSum[i] = 0;
        else rprefixSum[i] = tmpWeight[H - i - 1];
    }

    __syncthreads();

	for (int d = n / 2; d > 0; d /= 2)	// build sum in place up the tree 
	{ 	
	    if (id < d) { 
        	int ai = offset * (2 * id + 1) - 1;
        	int bi = offset * (2 * id + 2) - 1;
        
            if (rprefixSum[ai] != INT_MAX)
        	    rprefixSum[bi] += rprefixSum[ai];
        }
        
        offset *= 2;
	    __syncthreads();    
    }

    if (id == 0)
    {
    	rprefixSum[n - 1] = 0;
    }	// clear the last element  

    for (int d = 1; d < n; d *= 2)	// traverse down tree &build scan 
    { 	
        offset /= 2;      
        __syncthreads();      
        if (id < d) { 

        	int ai = offset * (2 * id + 1) - 1;
            int bi = offset * (2 * id + 2) - 1;
            
            int t = rprefixSum[ai];
            rprefixSum[ai] = rprefixSum[bi];
            if (t != INT_MAX)
                rprefixSum[bi] += t;
        }
    }
    __syncthreads();

    for (int i = id; i < H; i += stride) {
        prefixSum[i] = rprefixSum[H - i - 1];
    }

    __syncthreads();

    for (int i = id; i < H; i += stride){
        tmp[i] = sharedGridMap[i].x - prefixSum[i];
    }

    __syncthreads();

    // if (id == 0) {
    //     for (int i = 1; i < W; i++) {
    //         if (tmp[i] > tmp[i - 1]) {
    //             tmp[i] = tmp[i - 1];
    //             sharedGridMap[i].y = transTo1D(x, i - 1, W);
    //             *changed = 1;
    //         }
    //     }
    // }

    int startIdx = H - 1;
    for (int i = id; i < H; i += stride) {
        int ci = H - i - 1;
        if (ci >= 1) {
            for (int s = startIdx; s > ci; s--) {
                if (tmp[ci] > tmp[s]) {
                    tmp[ci] = tmp[s];
                    sharedGridMap[ci].y = transTo1D(ci + 1, y, W);
                    *changed = 1;
                }
            }
        }
        if (startIdx == H - 1) {
            startIdx -= (stride - 1);
        }
        else {
            startIdx -= stride;
        }
        __syncthreads();
    }

    __syncthreads();

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i].x = tmp[i] + prefixSum[i];
        gridMap[i * W + y] = sharedGridMap[i];
    }
}

__global__ void verticalSweep3(int H, int W, int2 *gridMap, int *hWeight, int *changed) {
    int y = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];

    if (id >= H || y >= W) return ;

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i] = gridMap[i * W + y];
        tmpWeight[i] = hWeight[i * W + y];
    }

    __syncthreads();

    for (int i = 1; i <= ceil(log2f(H)); i++) {
        for (int j = id; j < H; j += stride) {
            if (j % int(pow(2, i)) >= pow(2, i - 1)) {
                int idx = j - (j % int(pow(2, i - 1))) - 1;
                int partialCost = 0;

                for (int k = idx; k < j; k++) {
                    partialCost += tmpWeight[k];
                }

                if (sharedGridMap[idx].x != INT_MAX && sharedGridMap[j].x > sharedGridMap[idx].x + partialCost) {
                    sharedGridMap[j].x = sharedGridMap[idx].x + partialCost;
                    sharedGridMap[j].y = transTo1D(j - 1, y, W);
                    *changed = 1;
                }
            }
        }
    }

    __syncthreads();

    for (int i = id; i < H; i += stride) {
        gridMap[i * W + y] = sharedGridMap[i];
    }    
}

__global__ void verticalSweep4(int H, int W, int2 *gridMap, int *hWeight, int *changed) {
    int y = blockIdx.x;
    int id = threadIdx.x;
    int stride = blockDim.x;

    __shared__ int2 sharedGridMap[sharedSize];
    __shared__ int tmpWeight[sharedSize];

    if (id >= H || y >= W) return ;

    for (int i = id; i < H; i += stride) {
        sharedGridMap[i] = gridMap[i * W + y];
        tmpWeight[i] = hWeight[i * W  + y];
    }

    __syncthreads();

    for (int i = 1; i <= ceil(log2f(H)); i++) {
        for (int j = id; j < H; j += stride) {
            int invert_j = H - 1 - j;
            if (invert_j % int(pow(2, i)) >= pow(2, i - 1)) {
                int idx = H - 1 - (invert_j - (invert_j % int(pow(2, i - 1))) - 1);
                int partialCost = 0;

                for (int k = idx - 1; k >= j; k--) {
                    partialCost += tmpWeight[k];
                }

                if (sharedGridMap[idx].x != INT_MAX && sharedGridMap[j].x > sharedGridMap[idx].x + partialCost) {
                    sharedGridMap[j].x = sharedGridMap[idx].x + partialCost;
                    sharedGridMap[j].y = transTo1D(j + 1, y, W);
                    *changed = 1;
                }
            }
        }
    }

    __syncthreads();

    for (int i = id; i < H; i += stride) {
        gridMap[i * W + y] = sharedGridMap[i];
    }    
}

__global__ void getMinCostPin(int H, int W, int numOfPins, int2 *gridMap, int3 *pins) {
    int tid = threadIdx.x;
    int minPin = -1;
    // int minCost = INT_MAX;

    if (blockIdx.x != 0) return ;

    __shared__ int minPins[NumOfThreads];
    for (int i = threadIdx.x; i < NumOfThreads; i += blockDim.x) {
        if (i < numOfPins) minPins[i] = i;
        else minPins[i] = 0;
    }

    __syncthreads();

    for (int size = NumOfThreads / 2; size > 0; size /= 2) {
        if (tid < size) {
            int t = minPins[tid];
            int cmp = minPins[tid + size];
            if (pins[cmp].z == 0) {
                if (pins[t].z == 1 || gridMap[pins[cmp].x * W + pins[cmp].y].x < gridMap[pins[t].x * W + pins[t].y].x) {
                    minPins[tid] = cmp;
                }
            }
        }

        __syncthreads();
    }

    if (tid == 0) {
        minPin = minPins[0];
        // printf("%d\n", minPin);

        // retrace the path & set the cost = 0
        int x = pins[minPin].x;
        int y = pins[minPin].y;

        while (gridMap[x * W + y].x != 0) {
            // for debug
            if (gridMap[x * W + y].y == -1) {
                printf("getMinCost: cannot retrace the path\n");
                // cout << "getMinCost: cannot retrace the path" << endl;
            }

            gridMap[x * W + y].x = 0;
            int p = gridMap[x * W + y].y;
            x = p / W;
            y = p % W;
        }
        pins[minPin].z = 1;
    }


    // for testing 
    /*
    cout << "Pin: " << minPin << endl;
    cout << "Coord: " << pins[minPin]->x << " " << pins[minPin]->y << endl;
    cout << endl;
    */
}


// for testing
void printGridMap(int H, int W, int2 *gridMap) {
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (gridMap[i * W + j].x == INT_MAX) 
                cout << setw(3) << "INF" << " ";
            else 
                cout << setw(3) << gridMap[i * W + j].x << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// void sweep() {
//     int iteration = 0;
//     while (true) {
//         iteration++;
//         changed = false;
//         horizontalSweep();
//         verticalSweep();
        
//         if (!changed) break;
//     }

//     cout << "Iteration: " << iteration << endl;
// }

int main(int argc, char **argv) {
    // check the argument is enough
    if (argc != 3) {
        cout << "The number of arguments are wrong!!" << endl;
        return 0;
    }

    int W, H, numOfPins;
    int *hWeight, *vWeight;
    int2 *gridMap;
    int *changed;
    int3 *pins;

    globalTimer::startTimer();
    // parse basic
    parse_basic(argv[1], W, H, numOfPins);

    // new 
    hWeight = new int [W * H];
    vWeight = new int [W * H];
    gridMap = new int2 [W * H];
    changed = new int;
    pins = new int3 [numOfPins];

    // parse the input file
    parser(argv[1], W, H, numOfPins, pins, hWeight, vWeight, gridMap);
    std::chrono::nanoseconds _input = globalTimer::endTimer();

    // routing (sweep)
    gridMap[pins[0].x * W + pins[0].y].x = 0;
    pins[0].z = 1;

    // for testing 
    // printGridMap();

    // cuda
    const int numOfThreads = NumOfThreads;
    const int numOfBlocks = W;

    // cudaMalloc
    int *dhWeight, *dvWeight;
    int2 *dgridMap;
    int *dchanged;
    int3 *dpins;

    cudaMalloc(&dhWeight, H * W * sizeof(int));
    cudaMalloc(&dvWeight, H * W * sizeof(int));
    cudaMalloc(&dgridMap, H * W * sizeof(int2));
    cudaMalloc(&dpins, numOfPins * sizeof(int3));
    cudaMalloc(&dchanged, sizeof(int));

    cudaMemcpy(dhWeight, hWeight, H * W * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dvWeight, vWeight, H * W * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dgridMap, gridMap, H * W * sizeof(int2), cudaMemcpyHostToDevice);
    cudaMemcpy(dpins, pins, numOfPins * sizeof(int3), cudaMemcpyHostToDevice);
    cudaMemcpy(dchanged, changed, sizeof(int), cudaMemcpyHostToDevice);

    globalTimer::startTimer();
    for (int i = 0; i < numOfPins - 1; i++) {
        initialize<<<numOfBlocks, numOfThreads>>>(H, W, dgridMap);
        cudaDeviceSynchronize();

        int iteration = 0;
        while (true) {
            iteration++;
            *changed = 0;
            cudaMemcpy(dchanged, changed, sizeof(int), cudaMemcpyHostToDevice);

            horizontalSweep1<<<numOfBlocks, numOfThreads>>>(H, W, dgridMap, dvWeight, dchanged);
            cudaDeviceSynchronize();
            horizontalSweep2<<<numOfBlocks, numOfThreads>>>(H, W, dgridMap, dvWeight, dchanged);
            cudaDeviceSynchronize();
            verticalSweep3<<<numOfBlocks, numOfThreads>>>(H, W, dgridMap, dhWeight, dchanged);
            cudaDeviceSynchronize();
            verticalSweep4<<<numOfBlocks, numOfThreads>>>(H, W, dgridMap, dhWeight, dchanged);
            cudaDeviceSynchronize();

            cudaMemcpy(gridMap, dgridMap, H * W * sizeof(int2), cudaMemcpyDeviceToHost);
            printGridMap(H, W, gridMap);
            
            cudaMemcpy(changed, dchanged, sizeof(int), cudaMemcpyDeviceToHost);
            if (*changed == 0) break;
        }

        // cout << "Iteration: " << iteration << endl;

        // sweep();
        getMinCostPin<<<numOfBlocks, numOfThreads>>>(H, W, numOfPins, dgridMap, dpins);

        // for testing 
        // printGridMap();
    }
    std::chrono::nanoseconds _kernel = globalTimer::endTimer();

    globalTimer::startTimer();
    cudaMemcpy(gridMap, dgridMap, H * W * sizeof(int2), cudaMemcpyDeviceToHost);
    writeOutput(argv[2], H, W, gridMap);
    std::chrono::nanoseconds _output = globalTimer::endTimer();

    // cout time usage
    auto inputTime = _input.count();
    auto outputTime = _output.count();
    auto kernelTime = _kernel.count();

    // cout << "I/O time: " << inputTime / 1e9 + outputTime / 1e9 << endl;
    cout << "Kernel time: " << kernelTime / 1e9 << endl;

    // delete 
    delete []hWeight;
    delete []vWeight;
    delete []gridMap;
    delete []changed;
    delete []pins;

    // cudaFree
    cudaFree(dhWeight);
    cudaFree(dvWeight);
    cudaFree(dgridMap);
    cudaFree(dchanged);
    cudaFree(dpins);

    return 0;
}