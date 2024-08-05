#include "pivoter.h"
#include <cassert>

// #define DDEBUG
// #define BASELINE

#define BRANCHES

#ifdef BRANCHES
ui PNODES = 0, HNODES = 0;
#endif

#ifdef DDEBUG
ui uu = 0, uuu = 0;
#endif
std::vector<double> pivoter::run() {
    printf("pivoterSG.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
#endif
        ui deep = 1;
        nodes[deep].clear();
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            nodes[deep].push_back(g.pEdge[i]);
            // level[g.pEdge[i]] = deep;
        }

#ifdef DDEBUG
double p = answers[k];
#endif

if(nodes[deep].size() == 0) continue;
        
std::vector<ui> & C = nodes[deep];
std::vector<ui> & nxtC = nodes[deep + 1];
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
for(auto v : C) {
    if(!g.connectHash(maxV, v)) nxtC.push_back(v);
}
// std::vector<ui> cliqueNei, clique;
// clique.clear(); cliqueNei.clear();
// clique.push_back(maxV);
// for(auto v : C) {
//     if(g.connectHash(maxV, v)) cliqueNei.push_back(v);
//     else if(v != maxV) nxtC.push_back(v);
// }
// assert(C.size() == g.pIdx[u+1] - g.pIdx2[u]);
// if(nxtC.size() + clique.size() + cliqueNei.size() != C.size()) {
//     printf("%u %u %u %u\n", nxtC.size(), clique.size(), cliqueNei.size(), C.size());
//     fflush(stdout);
// }
// assert(nxtC.size() + clique.size() + cliqueNei.size() == C.size());

// while(cliqueNei.size() > 0) {
//     ui maxDeg = 0, maxV = cliqueNei[0];
//     findM(cliqueNei, maxDeg, maxV);
//     clique.push_back(maxV);

//     std::vector<ui> newCliqueNei;
//     for(ui i = 0; i < cliqueNei.size(); i++) {
//         if(g.connectHash(maxV, cliqueNei[i])) {
//             newCliqueNei.push_back(cliqueNei[i]);
//         }
//         else if(cliqueNei[i] != maxV) nxtC.push_back(cliqueNei[i]);
//     }
//     cliqueNei = newCliqueNei;
// }
// for(ui i = 0; i < clique.size(); i++) {
//     for(ui j = i + 1; j < clique.size(); j++) {
//         assert(g.connectHash(clique[i], clique[j]));
//     }
// }
// assert(nxtC.size() + clique.size() == C.size());
        deep += 1;
        for(ui v : nodes[deep]) sg.pIdx[v] = sg.deg[deep][v] = g.pIdx[v];
        for(ui v : nodes[deep]) level[v] = deep;
        for(ui v : nodes[deep]) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == deep) {
                    sg.pEdge[sg.deg[deep][v]++] = w;
                    sg.pEdge[sg.deg[deep][w]++] = v;
                }
            }
        }
        searchSG(deep, 0, 1);

        for(ui v : nodes[deep]) level[v] = 0;

#ifdef DDEBUG
printf("%u %.0f\n", u, answers[k] - p);
#endif
    }

#ifdef BASELINE
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != k) return false;

    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
        ui d = 0;
        for(ui v = 0; v < g.n; v++) if((1<<v) & x){
            if(g.connectHash(u, v)) d++;
        }
        if(d+1 < sz) return false; 
    }
    return true;
};

ui cnt = 0;
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        cnt++;
        print(i);
    }
}
printf("cnt:%u\n", cnt);

#endif

#ifdef BRANCHES
printf("hNodes:%u\npNodes:%u\n", HNODES, PNODES);
#endif
    return answers;
}

void pivoter::searchSG(ui deep, ui p, ui h) {
#ifdef DDEBUG
if(uu == uuu) {
printf("    deep %u\nC:", deep);
for(ui i = 0; i < nodes[deep].size(); i++) {
printf("%u ", nodes[deep][i]);
}printf("\n");
printf("p %u, h %u\n", p, h);
}
#endif
    if(h >= maxK) return;
    auto updateAns = [&]() {
        for(ui i = h; i < maxK && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
        }
    };
    std::vector<ui> & C = nodes[deep];
    ui edC = C.size();

    if(edC == 0) {
        updateAns(); return;
    }
    if(edC == 1) {
        p += 1; updateAns(); return;
    }
    if(edC == 2) {
        if(sg.deg[deep][C[0]] > sg.pIdx[C[0]]) {
            p += 2; updateAns();
        }
        else {
            p += 1; updateAns();
            p -= 1; h += 1; updateAns();
        }
        return;
    }

    ui pivot = C[0], pivotDeg = 0, num = 0;
    for(ui i = 0; i < edC; i++) {
        // ui tmp = 0;
        ui v = C[i];
        ui tmp = sg.deg[deep][v] - sg.pIdx[v];

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }
#ifdef DDEBUG
if(uu == uuu) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    if(pivotDeg+1 == edC && num == edC) {
        p += edC; updateAns(); return;
    }

    std::vector<ui> & nxtC = nodes[deep + 1];
    nxtC.clear();
    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui w = sg.pEdge[i];
        nxtC.push_back(w);
        level[w] = deep + 1;
    }
    auto updateSG = [&]() {
        for(ui j = 0; j < nxtC.size(); j++) {
            ui u = nxtC[j];
            ui & ed = sg.deg[deep + 1][u];
            ed = sg.deg[deep][u];
            for(ui l = sg.pIdx[u]; l < ed; ) {
                ui w = sg.pEdge[l];
                if(level[w] != deep + 1) {
                    std::swap(sg.pEdge[l], sg.pEdge[--ed]);
                }
                else l++;
            }
        }
    };
    updateSG();
#ifdef BRANCHES
PNODES++;
#endif
    searchSG(deep + 1, p + 1, h);
    for(auto v:nxtC) level[v] = deep;

#ifdef DDEBUG
if(uu == uuu) {
printf("    redeep %u, \nC:", deep);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
    if(pivotDeg+1 < C.size())
    for(ui u : C) {
        if(u == pivot) continue;
        if(g.connectHash(pivot, u)) continue;
        level[u] = deep - 1;
#ifdef DDEBUG
if(uu == uuu) {
printf("    Cdeep %u, cand %u\nC:", deep, u);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
        nxtC.clear();
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui w = sg.pEdge[j];
            if(level[w] == deep - 1) continue;
// assert(level[w] == deep);
            level[w] = deep + 1;
            nxtC.push_back(w);
        }
        updateSG();
#ifdef BRANCHES
HNODES++;
#endif
        searchSG(deep + 1, p, h + 1);

        for(auto v : nxtC) level[v] = deep;
    }
}

pivoter::pivoter(Graph && g, ui k) :g(g), k(k) { 
    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answers.resize(maxK, 0.0);

    computeC();

    sg.pIdx.resize(g.n);
    sg.pIdx2.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        sg.deg[i].resize(g.n);
    }

    nodes.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        nodes[i].resize(g.coreNumber);
    }

    level.resize(g.n);
    // ok.resize(g.n, false);

    printf("pivoter.h\n");
}

