#include "sdcCounting.h"
#include <queue>
#include <algorithm>
#include <utility>

// #define BASELINE

// #define DDEBUG

#include <iostream>

ui maxD = 0;
ui maxC = 0;
std::vector<double> sdcCounting::run() {
    printf("sdcCountingMembuffer.cpp::run");

    g.initHash();
    printf("init Hash\n");fflush(stdout);

    using Pair = std::pair<ui, ui>;
    std::queue<ui> que;

    for(ui u = 0; u < g.n; u++) {
        std::vector<ui> & C = nodes[0];
        C.clear();
// std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#ifdef DDEBUG
std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#endif
        //reduction to q-s-2 core
        std::vector<ui> & deg = neiInP;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];
            for(ui j = i + 1; j < g.pIdx[u + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connectHash(v, w)) deg[v]++, deg[w]++;
            }

            if(deg[v] < q-s-2) que.push(v);
        }
        while(!que.empty()) {
            ui v = que.front(); que.pop();

            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui w = g.pEdge[i];
                if(!g.connectHash(v, w)) continue;
                if((--deg[w]) == q-s-3) que.push(w);
            }
        }
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];
            if(deg[v] >= q-s-2) C.push_back(v);
            deg[v] = 0;
        }

        //2-hop
        ui edC1 = C.size();
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(w <= u || g.connectHash(u, w)) continue;

                if((++deg[w]) == q-s-2 || (q==s+2 && deg[w]==1)) {
                    C.push_back(w);
                }
            }
        }
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(w > u) deg[w] = 0;
            }
        }

        for(ui i = 0; i < edC1; i++) neiInP[C[i]] = 1;
        std::sort(C.begin(), C.end());
// for(auto v: C) printf("%u ", v);printf("\n");
        //build sub-graph g
        for(ui v : C) sg.pIdx[v] = sg.deg[0][v] = g.pIdx[v];
        for(ui v : C) level[v] = 1;
        level[u] = 1;
        auto buildSG = [&](ui v) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == 1) {
                    sg.pEdge[sg.deg[0][v]++] = w;
                    sg.pEdge[sg.deg[0][w]++] = v;
                }
            }
        };
        // for(ui v : C) buildSG(v);
        // buildSG(u);
        for(auto v : C) level[v] = 0;

        P.push_back(u);
        // listing(0, 0);
        P.pop_back();
        for(ui i = 0; i < C.size(); i++) neiInP[C[i]] = 0;
    }

#ifdef BASELINE
#define BASELINES 1
#define BASELINEQ 3
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto getMissEdges = [&](ui x) {
    ui TotalD = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        ui d = 0;
        for(ui v = u + 1 ; v < g.n; v++) if((1<<v) & x){
            if(!g.connectHash(u, v)) d++;
        }
        TotalD += d;
    }
    return TotalD;
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != BASELINEQ) return false;
    return getMissEdges(x) <= BASELINES;
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
    return answers;
}

void sdcCounting::pivoter(ui deep, ui p, ui h) {
// if(deep > maxD) {
//     maxD = deep;
//     printf("maxD:%u\n", maxD);
// }
    if(h == q) {
        answers[q] += 1;
        return;
    }

    std::vector<ui> & C = nodes[deep];
    auto updateAns = [&]() {
        for(ui i = h; i <= q && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
        }
    };

    if(C.size() == 0) {
        updateAns(); return;
    }
    if(C.size() == 1) {
        p += 1; updateAns(); return;
    }
    if(C.size() == 2) {
        if(g.connectHash(C[0], C[1])) {
            p += 2; updateAns();
        }
        else {
            p += 1; updateAns();
            p -= 1; h += 1; updateAns();
        }
        return;
    }

    ui pivot = C[0], pivotDeg = 0, pivotIndex = 0, num = 0;
    for(ui i = 0; i < C.size(); i++) {
        ui tmp = sg.deg[deep][C[i]] - sg.pIdx[C[i]];

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }

    if(pivotDeg+1 == C.size() && num == C.size()) {
        p += C.size(); updateAns(); return;
    }

    std::vector<ui> & nxtC = nodes[deep + 1];

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

    ui candSize = C.size() - pivotDeg;

    nxtC.clear();
    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui w = sg.pEdge[i];
        nxtC.push_back(w);
        level[w] = deep + 1;
    }

    updateSG();

    pivoter(deep + 1, p + 1, h);
    for(auto v:nxtC) level[v] = deep;

    if(pivotDeg+1 < C.size())
    for(ui u : C) {
        if(u == pivot) continue;
        if(g.connectHash(pivot, u)) continue;
        level[u] = deep - 1;

        nxtC.clear();
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui w = sg.pEdge[j];
            if(level[w] == deep - 1) continue;
            level[w] = deep + 1;
            nxtC.push_back(w);
        }
        updateSG();

        pivoter(deep + 1, p, h + 1);

        for(auto v : nxtC) level[v] = deep;
    }
}

void sdcCounting::listing(ui deep, ui missedEdges) {
    std::vector<ui> & C = nodes[deep];
#ifdef DDEBUG
printf("       deep %u\n", deep);
printf("C:");
for(ui i = 0; i < C.size(); i++) printf("%u ", C[i]); 
printf("neiInPC:");
for(ui i = 0; i < C.size(); i++) printf("%u ", neiInP[C[i]]); printf("\n");
printf("P:"); for(auto v:P) printf("%u ", v); printf("\n");
printf("missE %u\n", missedEdges);
printf("\n");
#endif
    
    if(P.size() == q - 1) {
#ifdef DDEBUG
printf("ans+%u\n", C.size());
#endif
        answers[q] += C.size();
        return;
    }
    if(missedEdges == s) {
        for(auto v : C) level[v] = deep;
        for(ui j = 0; j < C.size(); j++) {
            ui u = C[j];
            ui & ed = sg.deg[deep][u];
            ed = sg.deg[0][u];
            for(ui l = sg.pIdx[u]; l < ed; ) {
                ui w = sg.pEdge[l];
                if(level[w] != deep) {
                    std::swap(sg.pEdge[l], sg.pEdge[--ed]);
                }
                else l++;
            }
        }

        pivoter(deep, 0, P.size());

        for(auto v : C) level[v] = 0;
        return;
    }

    std::vector<ui> & nxtC = nodes[deep + 1];

    for(ui i = 0; i < C.size(); i++) {
        ui u = C[i];
        nxtC.clear();
#ifdef DDEBUG
printf("        deep %u, i %u, u %u\n", deep, i, u);
printf("C:");
for(ui i = 0; i < C.size(); i++) printf("%u ", C[i]); 
printf("neiInPC:");
for(ui i = 0; i < C.size(); i++) printf("%u ", neiInP[C[i]]); printf("\n");
printf("P:"); for(auto v:P) printf("%u ", v); printf("\n");
printf("missE %u\n", missedEdges);
#endif
        ui newMissEdges = missedEdges + P.size() - neiInP[u];
        // if(newMissEdges > s) continue;
        P.push_back(u);

        for(ui j = i + 1; j < C.size(); j++) {
            ui v = C[j];
            if(g.connectHash(u, v)) {
                if(newMissEdges + P.size() - neiInP[C[j]] <= s+1) {
                    neiInP[C[j]]++;
                    nxtC.push_back(C[j]);
                }
            }
            else {
                if(newMissEdges + P.size() - neiInP[C[j]] <= s) {
                    nxtC.push_back(C[j]);
                }
            }
        }
        
        listing(deep + 1, newMissEdges);

        for(auto v: nxtC) if(g.connectHash(u, v)) neiInP[v]--;
        P.pop_back();
    }
}


sdcCounting::sdcCounting(Graph && g, ui s, ui q):g(g), s(s), q(q) {
    ui maxDepth = (g.coreNumber + 1) * (s + 1) ;
    nodes.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        nodes[i].resize(g.n);
    }

    sg.pIdx.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        sg.deg[i].resize(g.n);
    }

    neiInP.resize(g.n);

    level.resize(g.n);

    answers.resize(q + 1);

    computeC();

    printf("sdcCounting.h\n");
}

sdcCounting::~sdcCounting() {
    delete [] bf3;
    delete [] CN;
    delete [] memBuffer;
}

