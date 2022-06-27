#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;

int main(int argc, char** argv) {
    // arguments check
    if (argc != 5) {
        cout << "The number of arguments are wrong!!" << endl;
        return 0;
    }

    int W = stoi(argv[1]);
    int H = stoi(argv[2]);
    int numOfPins = stoi(argv[3]);
    int weightRange = 90;

    ofstream fout(argv[4]);
    fout << W << " " << H << "\n";
    fout << numOfPins<< "\n";

    for (int i = 0; i < numOfPins; i++) {
        int x = rand() % H;
        int y = rand() % W;
        fout << "Pin " << x << " " << y << "\n";
    }

    for (int i = 0; i < H; i++) {
        fout << "Vertical ";
        for (int j = 0; j < W - 1; j++) {
            int weight = 1 + rand() % weightRange;
            fout << weight << " ";
        }
        fout << "\n";
    }

    for (int i = 0; i < H - 1; i++) {
        fout << "Horizontal ";
        for (int j = 0; j < W; j++) {
            int weight = 1 + rand() % weightRange;
            fout << weight << " ";
        }
        fout << "\n";
    }

    fout.close();
    return 0;
}