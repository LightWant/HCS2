#ifndef KCLISTCLIQUE_H
#define KCLISTCLIQUE_H

#include "../graph/graph.h"
#include "../tools/types.hpp"

class kclistClique {
private:
    Graph g;
    ui k = 0;
    double answer = 0;
    std::vector<std::vector<ui>> nodes;
    std::vector<std::vector<ui>> adj;
    std::vector<std::vector<ui>> edAdj;
    std::vector<uint8_t> level;
    std::vector<ui> clique;
    void listing(ui deep, ui edClique);

private:
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

public:
    kclistClique(Graph && g, ui k) :g(g), k(k) {
        maxSize = g.coreNumber + 1;
        maxK = k + 1;
        answer = 0;
        computeC();

        nodes.resize(maxK);
        for(ui i = 0; i < maxK; i++) {
            nodes[i].resize(maxSize);
        }

        edAdj.resize(k);
        for(ui i = 0; i < k; i++) edAdj[i].resize(g.n);
        adj.resize(g.n);
        for(ui i = 0; i < g.n; i++) {
            adj[i].resize(g.pIdx[i+1] - g.pIdx2[i]);
            for(ui j = 0; j < g.pIdx[i+1] - g.pIdx2[i]; j++) {
                adj[i][j] = g.pEdge[g.pIdx2[i] + j];
            }
        }

        // sg.pIdx.resize(g.n);
        // sg.pEdge.resize(g.m);
        // sg.deg.resize(maxK);
        // for(ui i = 0; i < maxK; i++) {
        //     sg.deg[i].resize(g.n);
        // }

        level.resize(g.n);

        clique.resize(g.coreNumber);
        clique.clear();

        printf("kCListClique.h\n");
    }
    ~kclistClique() { 
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
    }

    double run();
};

#endif