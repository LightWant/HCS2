#ifndef SDCLOCAL
#define SDCLOCAL

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <cassert>

class sdcLocal {
private:
    Graph g, sg;
    ui s = 0, q = 0;

    ui * candBuffer = nullptr;
    std::vector<ui> neiInH;
    std::vector<ui> level;
    std::vector<double> answers;

    std::vector<ui> P, rP;
    // std::vector<ui> cntNNeiInH;
    std::vector<std::vector<double>> dp;

private://edge idx
    std::vector<ui> pOEdge, pOIdx, reEegeId;
    std::vector<ui> Pw, H;
    std::vector<std::vector<ui>> binNodes;

private:
    void listing(ui deep, ui * C, ui sz, ui p, ui h, ui missEdges);

private://combinatorial number
    const ui maxSize = 1000;
    double ** CN = nullptr, *bf3 = nullptr;
    ui maxD = maxSize, maxD2 = maxSize;
    void computeC() {
        CN = new double*[maxD];
        bf3 = new double[maxD2 * maxD];
        for(int i = 0; i < maxD; i++) {
            CN[i] = bf3 + i * maxD2;
        }
        CN[0][0] = 1;
        CN[1][0] = 1;
        CN[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            CN[i][0] = 1;
            if(i < maxD2) CN[i][i] = 1;
            for(int j = 1; j < i && j < maxD2; j++) {
                CN[i][j] = CN[i - 1][j - 1] + CN[i - 1][j];
            }
        }
    }

private:
    ui * memBuffer = nullptr;
    ui * pointerBuffer = nullptr;
    void initBuffer(ui n) {
        memBuffer = new ui[n];
        pointerBuffer = memBuffer;
    }
    ui * allocMem(ui n) {
        ui * tmp = pointerBuffer;
        pointerBuffer += n;
        return tmp;
    }
    void freeMem(ui n) {
        pointerBuffer -= n;
    }

public:
    sdcLocal(Graph && g, ui s, ui q);
    ~sdcLocal() {
        if(bf3 != nullptr) delete [] bf3;
        if(CN != nullptr) delete [] CN;
        if(candBuffer != nullptr) delete [] candBuffer;
        if(memBuffer != nullptr) delete [] memBuffer;
    }

    void run();
};

#endif