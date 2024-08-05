#ifndef D2K_H
#define D2K_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include "candidate.hpp"
#include <vector>

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
    void pop(uint32_t u) { --cntNonNei[u]; }

    void clear(uint32_t u) { cntNonNei[u] = 0; }

    uint32_t * getBuffer(uint32_t u) {return buffer.data() + u * k;}
};

#define D2KNG
class d2k {
private:
    Graph g, sg;

    ui answer = 0;

    ui k = 0, q = 0;

    std::vector<ui> logBasedOn2;

    LinearSet C, X, P;

    nonNeiMatainer nn;//non-neighbors

private://the prune technique of dai
    std::vector<std::vector<ui>> bucket;
    std::vector<ui> support;

private:
    void bkPivot(ui deep, ui stX, ui edX, ui edC, ui edP);


#ifdef D2KNG
private://nadj
    std::vector<std::vector<ui>> nadj;
#endif

public:
    d2k(Graph && g, ui k, ui q):g(g), k(k), q(q) {  
        C.resize(g.n);
        X.resize(g.n);
        P.resize(g.n);

        sg.pIdx.resize(g.n);
        sg.pIdx2.resize(g.n);
        sg.pEdge.resize(g.m);

        nn.init(g.n, k);

        logBasedOn2.resize(g.maxD + 1);
        ui b = 0;
        for(ui i = 1; i <= g.maxD; i++) {
            while((1<<b) < i) b++;
            logBasedOn2[i] = b;
        }

        bucket.resize(k);
        support.resize(g.n);

        printf("dk2\n");
    }

    ui run();
};

#endif