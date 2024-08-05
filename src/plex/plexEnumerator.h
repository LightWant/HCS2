#ifndef PLEXENUMERATOR_H
#define PLEXENUMERATOR_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include "candidate.hpp"
#include <vector>

class plexEnumerator {
private:
    Graph g, sg;
    
    //key means the count of neighbors in P
    candidate neighborsInP;
    
    LinearSet P, C;
    ui edP, edC;
    
    LinearSet X;
    ui stX, edX;

    ui cnt = 0;
    ui k = 0, q = 0;

    //temp buffer
    std::vector<ui> nonNeighborsOfPv;

    void bkPivot(ui deep);

    std::vector<ui> logBasedOn2;

    std::vector<std::vector<ui>> bucket;
    std::vector<ui> support;

public:
    plexEnumerator(Graph && g, ui k, ui q):g(g), k(k), q(q) {  
        C.resize(g.n);
        edC = 0;
        P.resize(g.n);
        edP = 0;

        X.resize(g.n);
        stX = 0;
        edX = 0;

        neighborsInP.resize(g.n, g.maxD + 1);
        nonNeighborsOfPv.resize(g.n);

        sg.pIdx.resize(g.n);
        sg.pIdx2.resize(g.n);
        sg.pEdge.resize(g.m);

        cnt = 0;

        logBasedOn2.resize(g.maxD + 1);
        ui b = 0;
        for(ui i = 1; i <= g.maxD; i++) {
            while((1<<b) < i) b++;
            logBasedOn2[i] = b;
        }

        bucket.resize(k);
        support.resize(g.maxD);
    }

    ui run();
};

#endif