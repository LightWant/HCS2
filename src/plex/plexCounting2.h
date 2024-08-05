#ifndef PLEXCOUNTING2_H
#define PLEXCOUNTING2_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <cassert>

struct nonNeiMatainer{
    std::vector<uint32_t> buffer;
    std::vector<uint32_t> cntNonNei;
    uint32_t n;
    int k;

    void init(uint32_t n, int k) {
        this->n = n;
        this->k = k;
        buffer.resize(n * k);
        cntNonNei.resize(n, 0);
    }
    
    void addNonNei(uint32_t u, uint32_t nonNei) { 
        buffer[u * k + cntNonNei[u]++] = nonNei;
    }

    uint32_t getCntNonNei(uint32_t u) { return cntNonNei[u]; }
    void pop(uint32_t u) { 
// assert(cntNonNei[u] > 0);
        --cntNonNei[u]; 
    }

    void clear(uint32_t u) { cntNonNei[u] = 0; }

    uint32_t * getBuffer(uint32_t u) {return buffer.data() + u * k;}
};

class plexCounting {
private:
    Graph g, sg, nsg;
    ui k = 0, q = 0, Q = 0;
    LinearSet C, P;
    std::vector<std::vector<ui>> nodes;
    nonNeiMatainer nn;//non-neighbors
    std::vector<ull> answers;

private://combinatorial number
    const ui maxSize = 1000;
    double ** CN, *bf3;
    ui maxD = maxSize*10, maxD2 = maxSize;
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

private://nadj
    std::vector<std::vector<ui>> nadj;

private://the prune technique of dai
    std::vector<std::vector<ui>> bucket;
    std::vector<ui> support;

private:
    void bkPivot(ui deep, ui edC, ui edP, ui p, ui h);

public:
    plexCounting(Graph && g, ui k, ui q):g(g), k(k), q(q) {
        C.resize(g.n);
        P.resize(g.n);

        sg.pIdx.resize(g.n);
        sg.pEdge.resize(g.m);

        sg.deg.resize(g.n);
        sg.deg[0].resize(g.n);

        nn.init(g.n, k);

        bucket.resize(k);
        support.resize(g.n);

        answers.resize(maxSize, 0);

        computeC();

        printf("plexCounting\n");
    }
    ~plexCounting() {
        delete [] bf3;
        delete [] CN;
    }

    std::vector<ull> run();
};

#endif