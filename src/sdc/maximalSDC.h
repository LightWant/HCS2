#ifndef MAXIMALSDC_H
#define MAXIMALSDC_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <cassert>

class msdc {
private:
    Graph g, sg;
    ui s = 0, q = 0;
    LinearSet C, X;
    std::vector<ui> P;
    std::vector<ui> neiInP;
    ull answer = 0;

private:
    void bkPivot(ui deep, ui stX, ui edX, ui edC, ui missedEdges);

public:
    msdc(Graph && g, ui s, ui q):g(g), s(s), q(q) {
        C.resize(g.n);
        X.resize(g.n);

        sg.pIdx.resize(g.n);
        sg.pIdx2.resize(g.n);
        sg.pEdge.resize(g.m);

        sg.deg.resize(g.n);

        neiInP.resize(g.n);

        // bucket.resize(s);
        // support.resize(g.n);

        printf("msdc.h\n");
    }
    ull run();
};

#endif