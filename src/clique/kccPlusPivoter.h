#ifndef KCCPLUSPIVOTER_H
#define KCCPLUSPIVOTER_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"

class kccPlusPivoter {
private:
    Graph g;
    ui k = 0;
    LinearSet C;
    double answer = 0;

private://kclist
    std::vector<std::vector<ui>> nodes;
    std::vector<ui> level;

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

private:
    Graph sg;
    std::vector<ui> clique, Ps, Hs;
    void searchSGClique(ui deep, ui edC, ui p, ui h);
    
    std::vector<std::vector<ui>> L, R;
    std::vector<ui> weight;
    std::vector<ui> twoHopNodes;
    std::vector<bool> vis;
    std::vector<bool> inP;

    std::vector<std::vector<ui>> degL, degR;

private:
    void listClique(ui deep, ui edClique);

public:
    kccPlusPivoter(Graph && g, ui k);
    ~kccPlusPivoter() { 
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
    }

    double run();

};

#endif