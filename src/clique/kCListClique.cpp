#include "kCListClique.h"
#include <algorithm>

// #define DDEBUG

double kclistClique::run() {
    printf("kClistClique.cpp\n");
    g.initHash();
    printf("initHash\n");

#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] > 0) {
        nodes[0].clear();
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            nodes[0].push_back(g.pEdge[i]);
            level[g.pEdge[i]] = 1;
        }
#ifdef DDEBUG
printf("node[0]:");
for(auto v:nodes[0]) printf("%u ", v); printf("\n");
#endif
        std::vector<ui> & C = nodes[0];
        std::vector<ui> & nxtC = nodes[1];
        nxtC.clear();

        ui maxDeg = 0, maxV = C[0];
        auto findM = [&](std::vector<ui> & C, ui & maxDeg, ui & maxV) {
            for(ui i = 0; i < C.size(); i++) {
                ui deg = 0;
                for(ui j = 0; j < C.size(); j++) {
                    if(g.connectHash(C[i], C[j])) deg++;
                }
                if(deg > maxDeg) {
                    maxDeg = deg;
                    maxV = C[i];
                }
            }
        };
        findM(C, maxDeg, maxV);
        std::vector<ui> cliqueNei;
        clique.clear(); cliqueNei.clear();
        clique.push_back(maxV);
        for(auto v : C) {
            if(g.connectHash(maxV, v)) cliqueNei.push_back(v);
            else if(v != maxV) nxtC.push_back(v);
        }

        while(cliqueNei.size() > 0) {
            ui maxDeg = 0, maxV = cliqueNei[0];
            findM(cliqueNei, maxDeg, maxV);
            clique.push_back(maxV);

            std::vector<ui> newCliqueNei;
            for(ui i = 0; i < cliqueNei.size(); i++) {
                if(g.connectHash(maxV, cliqueNei[i])) {
                    newCliqueNei.push_back(cliqueNei[i]);
                }
                else if(cliqueNei[i] != maxV) nxtC.push_back(cliqueNei[i]);
            }
            cliqueNei = newCliqueNei;
        }

        std::sort(nxtC.begin(), nxtC.end());
        for(auto u : clique) level[u] = 0;
#ifdef DDEBUG
for(auto v:nxtC) assert(level[v] == 1);
assert(clique.size() + nxtC.size() == C.size());
for(ui i = 0; i < g.n; i++) printf("level %u %u\n", i, level[i]);
#endif
        for(auto v:nxtC) {
            ui ed = adj[v].size();
            for(ui i = 0; i < ed; ) {
                if(level[adj[v][i]] != 1) std::swap(adj[v][i], adj[v][--ed]);
                else i++;
            }
            edAdj[1][v] = ed;
#ifdef DDEBUG
for(ui i = 0; i < edAdj[1][v]; i++) assert(level[adj[v][i]] == 1);
#endif
        }
#ifdef DDEBUG
printf("subgraph %u\n", u);
for(auto v:nxtC) {
    printf("%u:", v);
    for(ui i = 0; i < edAdj[1][v]; i++) printf("%u ", adj[v][i]);
    printf("\n");
}
#endif
        listing(1, clique.size());

        for(ui v : nxtC) level[v] = 0;
    }

    return answer;
}

void kclistClique::listing(ui deep, ui edClique) {
#ifdef DDEBUG
printf("    deep %u\n", deep);
printf("C:");
for(auto v: nodes[deep]) printf("%u ", v);printf("\n");
printf("clique:");
for(ui i = 0; i < edClique; i++) printf("%u ", clique[i]);printf("\n");
#endif

    std::vector<ui> & C = nodes[deep];
    if(deep == k-1) {
        // for(auto v:C) answer += edAdj[deep][v];;
        answer += C.size();
        answer += edClique;

        return;
    }

    if(edClique >= k - deep)
        answer += CN[edClique][k - deep];

    std::vector<ui> & nxtC = nodes[deep + 1];
    
    for(ui i = 0; i < C.size(); i++) {
        ui u = C[i];
        if(edAdj[deep][u]+deep+edClique+1 < k) continue;

        nxtC.clear();
        for(ui j = 0; j < edAdj[deep][u]; j++) {
            ui v = adj[u][j];
            nxtC.push_back(v);
            level[v] = deep + 1;
        }

        for(ui j = 0; j < nxtC.size(); j++) {
            ui v = nxtC[j];
            ui & ed = edAdj[deep + 1][v];
            ed = edAdj[deep][v];

            for(ui l = 0; l < ed; ) {
                ui w = adj[v][l];
                if(level[w] == deep + 1) l++;
                else std::swap(adj[v][l], adj[v][--ed]);
            }
        }

        ui newEdClique = edClique;
        for(ui j = 0; j < newEdClique; ) {
            if(!g.connectHash(u, clique[j])) {
                std::swap(clique[j], clique[--newEdClique]);
            }
            else j++;
        }

        listing(deep + 1, newEdClique);

        for(auto v : nxtC) level[v] = deep;
    }
}
