#ifndef SDCCOUNTING2_H
#define SDCCOUNTING2_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <cassert>

class sdcCounting
{
private:
    Graph g, sg;
    ui s = 0, q = 0, Q = 0;
    
    ui * candBuffer = nullptr;
    std::vector<ui> neiInPH;
    std::vector<ui> level;
    std::vector<double> answers;

private:
    void listing(ui deep, ui * C, ui cm, ui sz, ui p, ui h, ui missEdges);
    void pivoter(ui deep, ui * C, ui sz, ui p, ui h);

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
    sdcCounting(Graph && g, ui s, ui q, ui Q);
    ~sdcCounting();

    std::vector<double>  run();
};




#endif