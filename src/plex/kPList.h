#ifndef KPLIST_H
#define KPLIST_H

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

class kplist {
private:
    Graph g, sg;
    ui k = 0, q = 0;
    
    //std::vector<ui> P;
    LinearSet P;
    ui edP = 0;
    std::vector<ui> level;
    nonNeiMatainer nn;//non-neighbors
    ull answer = 0;

private://the prune technique of dai
    std::vector<std::vector<ui>> bucket;
    std::vector<ui> support;

private:
    void listing(ui deep, const std::vector<ui> & C);

public:
    kplist(Graph && g, ui k, ui q):g(g), k(k), q(q) { 
        // C.resize(g.n);
        P.resize(g.n);

        sg.pIdx.resize(g.n);
        sg.pIdx2.resize(g.n);
        sg.pEdge.resize(g.m);

        nn.init(g.n, k);

        bucket.resize(k);
        support.resize(g.n);

        level.resize(g.n, 0);

        printf("kPlist\n");fflush(stdout);
    }

    ull run();
};

#endif