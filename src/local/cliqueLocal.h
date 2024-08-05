#ifndef CLIQUELOCAL_H
#define CLIQUELOCAL_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <algorithm>
#include <set>
#include <cassert>

class cliqueLocal {
private:
    Graph g;
    ui k = 0;
    LinearSet C;
    std::vector<double> answers;

private://edge idx
    std::vector<ui> pOEdge, pOIdx, reEegeId;

private://combinatorial number
    ui maxK = 31;
    ui maxSize = 3000;
    double ** CN = nullptr, *bf3 = nullptr;
    void computeC() {
        ui maxD = maxSize, maxD2 = maxK;
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

private://raw pivoter
    std::vector<ui> P, H, PP, PH, HH;
    void search(ui edC, ui p, ui h);
    
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
    cliqueLocal(Graph && g, ui k);
    ~cliqueLocal() { 
        if(bf3 != nullptr) {
            delete [] bf3; 
            delete [] CN; 
            bf3 = nullptr;
        }
        if(memBuffer != nullptr) {
            delete [] memBuffer; 
            memBuffer = nullptr;
        }
    }

    void run();

};

#endif