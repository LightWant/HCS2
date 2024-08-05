#include "sdcCounting2DP.h"
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
    printf("sdcCounting2DPCnt.cpp::run");

    g.initHash();
    printf("init Hash\n");fflush(stdout);

    using Pair = std::pair<ui, ui>;
    std::vector<ui> que(g.coreNumber * g.coreNumber);
    std::vector<ui> deg(g.n);
    ui lq, rq;
    std::vector<bool> vis(g.n);
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
        if(g.pIdx[u + 1] - g.pIdx2[u] + s + 1 < q) continue;

        ui * C = candBuffer;
        ui sz = 0;
#ifdef DDEBUG
uu = u;
#endif

#ifdef DDEBUG
std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#endif
        //reduction to q-s-2 core
        lq = rq = 0;
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
            if(deg[v] >= q-s-2) C[sz++] = v;
            deg[v] = 0;
        }
#ifdef DDEBUG
printf("edC1 %u:", sz);
for(ui i = 0; i < sz; i++) printf("%u ", C[i]); printf("\n");
#endif

        //2-hop
        ui edC1 = sz;
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            if(g.pIdx[v+1] > g.pIdx[v])
            for(ui j = g.pIdx[v+1] - 1; j >= g.pIdx[v]; j--) {
                ui w = g.pEdge[j];

                if(w == u) break;
                if(vis[w]) continue;

                if(++deg[w] == q-s-1) {//non-nei have at least q-s-1 common neighbors
                    C[sz++] = w;
                    vis[w] = true;
                }

                if(j == 0) break;
            }
        }
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];

            if(g.pIdx[v+1] > g.pIdx[v])
            for(ui j = g.pIdx[v+1] - 1; j >= g.pIdx[v]; j--) {
                ui w = g.pEdge[j];
                if(w == u) break;

                deg[w] = 0;
                if(j == 0) break;
            }
        }
        for(ui i = edC1; i < sz; i++) {
            vis[C[i]] = false;
        }
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            vis[g.pEdge[i]] = false;
        }

        for(ui i = 0; i < edC1; i++) neiInH[C[i]] = 1;
        std::sort(C + edC1, C + sz);

        //build sub-graph g
        for(ui i = 0; i < sz; i++) {
            ui v = C[i];
            sg.pIdx[v] = sg.deg[1][v] = g.pIdx[v];
        }
        for(ui i = 0; i < sz; i++) level[C[i]] = 1;
        // level[u] = 1;
        auto buildSG = [&](ui v) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == 1) {
                    sg.pEdge[sg.deg[1][v]++] = w;
                    sg.pEdge[sg.deg[1][w]++] = v;
                }
            }
        };
        for(ui i = 0; i < sz; i++) buildSG(C[i]);
        // buildSG(u);

#ifdef DDEBUG
double preAns = answers[Q];
#endif

        P.clear();
        listing(1, C, sz, 1, 0);

        for(ui i = 0; i < sz; i++) neiInH[C[i]] = 0;
        for(ui i = 0; i < sz; i++) level[C[i]] = 0;
        level[u] = 0;

#ifdef DDEBUG
printf("diff %u %.0f\n", u, answers[Q] - preAns);
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
    if(sz != Q) return false;
    return getMissEdges(x) <= s;
};

ui cnt = 0;
std::vector<ui> cnti(g.n);
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        cnt++;
        print(i);
        for(ui u = 0; u < g.n; u++) if((1<<u) & i) {
            cnti[u]++; break;
        }
    }
}
printf("cnt:%u\n", cnt);
for(ui u = 0; u < g.n; u++) {
    printf("cnt%u:%u\n", u, cnti[u]);
}
#endif
    return answers;
}

//C, cm, sz, h, missEdges, neiInPH, level, P
void sdcCounting::listing(ui deep, ui * C, ui sz, ui h, ui missEdges) {
#ifdef DDEBUG
printf("       deep %u, h %u, missE %u\n", deep, h, missEdges);
printf("C_sz:");
for(ui i = 0; i < sz; i++) printf("%u ", C[i]); 
printf("neiInH:");
for(ui i = 0; i < sz; i++) printf("%u ", neiInH[C[i]]); printf("\n");
printf("cntNeiInH:"); 
for(ui i = 0; i <= s; i++) printf("%u ", cntNNeiInH[i]); printf("\n");

printf("\n");
#endif

    ui rp = 0;
    for(ui i = 1; i <= s-missEdges; i++) rp += cntNNeiInH[i];
    if(h == Q-1) {
// #ifdef DDEBUG
// printf("ans+%u\n", sz+p);
// #endif
        answers[Q-1] += 1;
        answers[Q] += sz + rp + cntNNeiInH[0];
#ifdef DDEBUG
printf("ans0 + %u %u %u\n", sz, p, rp);
#endif
        return;
    }
    if(h+sz+rp+cntNNeiInH[0] < q) return;

    auto updateAns = [&]() {
#ifdef DDEBUG
printf("updateAns, h %u\n", h);
printf("cntNNeiInH:"); 
for(ui i = 0; i <= s; i++) printf("%u ", cntNNeiInH[i]); printf("\n");
#endif
        
        for(ui j = 1; j <= rp; j++) {
            for(ui k = 0; k <= s-missEdges; k++) {
                dp[j][k] = 0;
            }
        }
        for(ui i = 1; i <= rp; i++) {
            for(ui j = i; j >= 1; j--) {
                for(ui k = rP[i-1]; k <= s-missEdges; k++) {
                    dp[j][k] += dp[j-1][ k-rP[i-1] ];
                }
            }
        }

        ui p = cntNNeiInH[0];
        for(ui j = 0; j <= rp; j++) {
            for(ui k = 0; k <= s-missEdges; k++) {
                for(ui l = 0; l <= p; l++) {
                    if(h+j+l > Q) break;
#ifdef DDEBUG
if(h+j+l == Q) printf("ans + %.0f\n", dp[j][k]*CN[p][l]);
#endif
                    answers[h+j+l] += dp[j][k]*CN[p][l];
                }
            }
        }
    };
    if(sz == 0) {
        updateAns();
        
        return;
    }
    if(sz == 1) {
        cntNeiInH[ h-neiInH[C[0]] ]++;
        
        updateAns();
        
        cntNeiInH[ h-neiInH[C[0]] ]--;

        return;
    }

    ui * nxtC = C + sz;
    auto updateSG = [&](ui newSz) {
        for(ui j = 0; j < newSz; j++) {
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
    auto addPNode = [&](ui u, ui & newSz) {
        for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
            ui v = sg.pEdge[i]; 
            // neiInH[v]++;
            level[v] = deep+1;
            nxtC[newSz++] = v;
        }

        updateSG(newSz);
    };

    auto addHNode = [&](ui u, ui & newSz, ui newMissEdges) {
        for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
            ui v = sg.pEdge[i]; 
            if(level[v] != deep) continue;
            neiInH[v]++;
        }

        for(ui i = 0; i < sz; i++) {
            ui v = C[i];
            if(level[v] < deep) continue;

            if(newMissEdges + h+1-neiInH[v] <= s) {
                level[v] = deep+1;
                nxtC[newSz++] = v;
            }
        }

        updateSG(newSz);
    };


    ui pivot = C[0], pivotDeg = 0;
    for(ui i = 0; i < sz; i++) {
        ui u = C[i];
        ui deg = sg.deg[deep][u] - sg.pIdx[u];

        if(deg > pivotDeg) {
            pivotDeg = deg;
            pivot = u;
        }
    }
#ifdef DDEBUG
printf("pivot %u, deep %u\n", pivot, deep);
#endif
    //pivot branch
    ui newSz = 0;
    addPNode(pivot, newSz);

    ui candSise = sz - pivotDeg - 1;
    ui * cand = allocMem(candSise);
    for(ui i = 0, j = 0; i < sz; i++) {
        ui v = C[i];
        if(level[v] != deep + 1 && v != pivot) {
            cand[j++] = v;
// assert(j <= candSise);
        }
    }

    if(neiInH[pivot] == h)
        listing(deep+1, nxtC, newSz, p+1, h, missEdges);
    else {
        P.push_back(h-neiInH[pivot]);
        listing(deep+1, nxtC, newSz, p, h, missEdges);
        P.pop_back();
    }

    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui v = sg.pEdge[i]; 
        // neiInH[v]--;
        level[v] = deep;
    }
#ifdef DDEBUG
printf("       After deep %u, p %u, h %u, missE %u\n", deep, p, h, missEdges);

printf("C_sz:");
for(ui i = 0; i < sz; i++) printf("%u ", C[i]); 
printf("neiInH:");
for(ui i = 0; i < sz; i++) printf("%u ", neiInH[C[i]]); printf("\n");
printf("cands:");
for(ui i = 0; i < candSise; i++) printf("%u ", cand[i]); printf("\n");

printf("\n");
#endif

    for(ui i = 0; i < candSise; i++) {
        ui u = cand[i];
#ifdef DDEBUG
printf("  deep %u, cand %u, neiIPH %u\n", deep, u, neiInH[u]);
printf("cands:");
for(ui i = 0; i < candSise; i++) printf("%u ", cand[i]); printf("\n");
#endif
        level[u] = deep-1;
        ui newSz = 0, newMissEdges = missEdges + h-neiInH[u];
        addHNode(u, newSz, newMissEdges);

        listing(deep+1, nxtC, newSz, p, h+1, newMissEdges);

        for(ui j = 0; j < newSz; j++) level[nxtC[j]] = deep;
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui v = sg.pEdge[j]; 
            if(level[v] != deep) continue;
            neiInH[v]--;
        }
    }

    freeMem(candSise);
}

sdcCounting::sdcCounting(Graph && g, ui s, ui q, ui Q):g(g), s(s), q(q), Q(Q) {
    ui maxDepth = g.coreNumber + s + 5;
    ui maxCSize = std::min(g.coreNumber * g.coreNumber, g.n);
    candBuffer = new ui[maxDepth * maxCSize];

    initBuffer(maxDepth * maxCSize);

    sg.pIdx.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        sg.deg[i].resize(g.n);
    }

    neiInH.resize(g.n);

    level.resize(g.n);

    answers.resize(Q + 1);

    dp.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        dp[i].resize(s + 1);
    }
    dp[0][0] = 1;

    cntNNeiInH.resize(s + 1);

    maxD = maxDepth;
    maxD2 = Q + 1;
    computeC();

    printf("sdcCounting2DPCnt::sdcCounting\n");
}

sdcCounting::~sdcCounting() {
    if(bf3 != nullptr) delete [] bf3;
    if(CN != nullptr) delete [] CN;
    if(candBuffer != nullptr) delete [] candBuffer;
    if(memBuffer != nullptr) delete [] memBuffer;
}

