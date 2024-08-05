#ifndef D2KRAW_H
#define D2KRAW_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include "candidate.hpp"
#include <vector>

class d2kRaw {
private:
    Graph g, sg, sng;
    ui answer = 0;
    ui k = 0, q = 0;
    
    std::vector<ui> nn;
    void addNN(std::vector<ui> & cnn, ui w, ui u) {
        nn[w*k + cnn[w]++] = u;
    }

private://the prune technique of dai
    std::vector<std::vector<ui>> bucket;
    std::vector<ui> support;

private:
    using vc = std::vector<ui>;
    void bkPivot(ui deep, vc & P, vc & C, vc & X, vc & cnn);

public:
    d2kRaw(Graph && g, ui k, ui q):g(g), k(k), q(q) {  
        sg.pIdx.resize(g.n);
        sg.pIdx2.resize(g.n);
        sg.pEdge.resize(g.m);

        nn.resize(g.n * k);

        bucket.resize(k);
        support.resize(g.n);

        printf("dk2_raw\n");
    }

    ui run();

};

#endif