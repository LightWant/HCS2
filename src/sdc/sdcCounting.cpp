#include "sdcCounting.h"
#include <queue>
#include <algorithm>
#include <utility>

// #define BASELINE

// #define DDEBUG

#include <iostream>

#ifdef DDEBUG
ui uu = 0, uuu = 0;
#endif

ui maxD = 0;
ui maxC = 0;
std::vector<double> sdcCounting::run() {
    printf("sdcCounting.cpp::run");

    g.initHash();
    printf("init Hash\n");fflush(stdout);

    using Pair = std::pair<ui, ui>;
    std::vector<ui> que(g.coreNumber * g.coreNumber);
    ui lq, rq;
    std::vector<bool> vis(g.n);
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
        std::vector<ui> & C = nodes[0];
        C.clear();
#ifdef DDEBUG
uu = u;
#endif
// std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#ifdef DDEBUG
std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#endif
        //reduction to q-s-2 core
        std::vector<ui> & deg = neiInP;
        lq = rq = 0;
// for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
//     vis[g.pEdge[i]] = true;
// }
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) deg[g.pEdge[i]] = 0;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];
            for(ui j = i + 1; j < g.pIdx[u + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connectHash(v, w)) deg[v]++, deg[w]++;
            }

            if(deg[v] < q-s-2) {
                que[rq++] = v;
            }
        }
        while(lq < rq) {
            ui v = que[lq++];

            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui w = g.pEdge[i];
                if(!g.connectHash(v, w)) continue;
                if((--deg[w]) == q-s-3) que[rq++] = w;
            }
        }
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];

            vis[v] = true;
            if(deg[v] >= q-s-2) {
                C.push_back(v);
            }
            deg[v] = 0;
        }
#ifdef DDEBUG
printf("edC1 %u:", C.size());
for(auto v : C) printf("%u ", v); printf("\n");
#endif
// assert(C.size() <= g.pIdx[u+1] - g.pIdx2[u]);
        //2-hop
        ui edC1 = C.size();
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            if(g.pIdx[v+1] > g.pIdx[v])
            for(ui j = g.pIdx[v+1] - 1; j >= g.pIdx[v]; j--) {
                ui w = g.pEdge[j];
// printf("check %u-%u\n", v, w);
                if(w == u) break;
                if(vis[w]) continue;

                // if((++deg[w]) == q-s-2 || (q==s+2 && deg[w]==1)) {
                if(++deg[w] == q-s-1) {//non-nei have at least q-s-1 common neighbors
                    C.push_back(w);
                    vis[w] = true;
                }

                if(j == 0) break;
            }
        }
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
// printf("uncheck %u, %u %u\n", v, g.pIdx[v+1], g.pIdx[v]);
            if(g.pIdx[v+1] > g.pIdx[v])
            for(ui j = g.pIdx[v+1] - 1; j >= g.pIdx[v]; j--) {
                ui w = g.pEdge[j];
// printf("uncheck %u-%u\n", v, w);
                if(w == u) break;

                deg[w] = 0;
                if(j == 0) break;
            }
        }
        for(ui i = edC1; i < C.size(); i++) {
            vis[C[i]] = false;
        }
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            vis[g.pEdge[i]] = false;
        }
#ifdef DDEBUG
for(auto v:C) if(neiInP[v] != 0) printf("no 0: %u\n", v);
// for(auto v:C) assert(neiInP[v] == 0);
#endif 
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
        for(ui v : C) buildSG(v);
        buildSG(u);
        for(auto v : C) level[v] = 0;
        level[u] = 0;

// for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
//     vis[g.pEdge[i]] = true;
// }
#ifdef DDEBUG
double preAns = answers[q];
#endif
        P.push_back(u);
        listing(0, 0);
        P.pop_back();
        for(ui i = 0; i < C.size(); i++) neiInP[C[i]] = 0;

#ifdef DDEBUG
printf("diff %u %.0f\n", u, answers[q] - preAns);
#endif
    }

#ifdef BASELINE
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto getMissEdges = [&](ui x) {
    ui TotalD = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        ui d = 0;
        for(ui v = u + 1 ; v < g.n; v++) if((1<<v) & x) {
// assert(u < g.n && v < g.n);
// printf("find %u %u\n", u, v);
            if(!g.connectHash(u, v)) d++;
        }
        TotalD += d;
    }
    return TotalD;
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != q) return false;
    return getMissEdges(x) <= s;
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
#ifdef DDEBUG
printf("pdeep %u, p %u, h %u\nC:", deep, p, h);
for(auto v : nodes[deep]) printf("%u ", v); printf("\n");
#endif

    if(h == q-1) {
        answers[q] += p + nodes[deep].size();
#ifdef DDEBUG
printf("ans += %u, pivoter1\n", p + nodes[deep].size());
#endif
        return;
    }

    std::vector<ui> & C = nodes[deep];
    auto updateAns = [&]() {
        for(ui i = h; i <= q && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
#ifdef DDEBUG
if(i==q)
printf("ans += C[%u][%u] %.0f, pivoter2\n", p, i-h, CN[p][i-h]);
#endif
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
#ifdef DDEBUG
printf("pivot %u, pivotdeg %u, num %u\n", pivot, pivotDeg, num);
#endif
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
printf("neiInP:");
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
    // if(missedEdges == s &&  && C.size() > 8) {
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
#ifdef DDEBUG
printf("pivoter parameter h %u\n", P.size());
printf("subgrah:");
for(auto v: C) {
    printf("%u:", v);
    for(ui i = sg.pIdx[v]; i < sg.deg[deep][v]; i++) {
        ui w = sg.pEdge[i];
        printf("%u ", w);
    }printf("\n");
}
#endif
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
}

