#include "sdcLocal.h"
#include <queue>
#include <algorithm>
#include <utility>

// #define BASELINE
// #define DDEBUG

#ifdef DDEBUG
#include <iostream>
#endif

// #define DEBUG2

#ifdef DEBUG2
ui uu = 0;
ui uuu = 2089580;
#endif
void sdcLocal::run() {
    printf("sdcVLocal.cpp::run\n");fflush(stdout);

    g.initHash();
    // printf("init Hash\n");fflush(stdout);

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
#ifdef DEBUG2
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

        H.push_back(u);
        listing(1, C, sz, 0, 1, 0);
        H.pop_back();

        for(ui i = 0; i < sz; i++) neiInH[C[i]] = 0;
        for(ui i = 0; i < sz; i++) level[C[i]] = 0;
        level[u] = 0;

    }

    // for(ui i = 0; i < g.m/2; i++) {
    //     printf("%.0f\n", answers[reEegeId[i]]);
    // }

#ifdef BASELINE
std::vector<ui> ec(g.m/2, 0);
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");

    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        for(ui v = u + 1; v < g.n; v++) if((1<<v) & x) {
            if(!g.connectHash(u, v)) continue;
            auto st = pOEdge.begin() + pOIdx[u];
            auto ed = pOEdge.begin() + pOIdx[u + 1];
            ui idx = std::lower_bound(st, ed, v) - pOEdge.begin(); 
            ec[idx]++;
        }
    }
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
for(ui e = 0; e < g.m/2; e++) {
    printf("e_cnt %u, %u %.0f\n", e, ec[e], answers[e]);
}
for(ui u = 0; u < g.n; u++) {
    for(ui e = pOIdx[u]; e < pOIdx[u+1]; e++) {
        printf("e%u, %u %u\n", e, u, pOEdge[e]);
    }
}
#endif
}

void sdcLocal::listing(ui deep, ui * C, ui sz, ui p, ui h, ui missEdges) {

#ifdef DDEBUG
printf("       deep %u, p %u, h %u, missE %u\n", deep, p, h, missEdges);
printf("C_sz:");
for(ui i = 0; i < sz; i++) printf("%u ", C[i]); 
printf("neiInH:");
for(ui i = 0; i < sz; i++) printf("%u ", neiInH[C[i]]); printf("\n");
printf("H:"); for(auto v : H) printf("%u ", v); printf("\n");
for(ui i = 0; i <= s; i++) {
    printf("bin%u:", i);
    for(auto v : binNodes[i]) printf("%u ", v);
    printf("\n");
}
printf("\n");
#endif
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1\n");fflush(stdout);
}
#endif
    if(h == q-1) {
        for(ui i = 0; i < sz; i++) {
            answers[C[i]]++;
#ifdef DDEBUG
printf("vc %u, add 1, h==q-1\n", C[i]);
#endif  
        }
        double hAddV = sz;
        for(ui miss = 0; miss <= s-missEdges; miss++) {
            for(auto u : binNodes[miss]) {
                answers[u]++;
            }
            hAddV += binNodes[miss].size();
        }
        
        for(ui i = 0; i < H.size(); i++) {
            answers[C[i]] += hAddV;
#ifdef DDEBUG
printf("vh %u, add 1, h==q-1\n", C[i]);
#endif  
        }

        return;
    }
    if(h+p+sz < q) return;

    auto updateAns = [&]() {
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.2\n");fflush(stdout);
}
#endif
#ifdef DDEBUG
printf("updateAns\n");
for(ui i = 0; i <= s; i++) {
    printf("bin%u:", i);
    for(auto v : binNodes[i]) printf("%u ", v);
    printf("\n");
}
printf("H:");
for(auto v : H) printf("%u ", v); printf("\n");
#endif  
        //addH
        ui sumNodes = 0;
        for(ui i = 0; i <= s-missEdges; i++) sumNodes += binNodes[i].size();
        for(ui i = 0; i <= sumNodes; i++) {
            for(ui k = 0; k <= s-missEdges; k++) {
                dp[i][k] = 0;
            }
        }
        for(ui miss = 1, i = 0; miss <= s-missEdges; miss++) {
            for(ui x = 0; x < binNodes[miss].size(); x++) {
                i++;
                for(ui j = i; j >= 1; j--) {
                    for(ui k = miss; k <= s-missEdges; k++) {
                        dp[j][k] += dp[j-1][ k-miss ];
                    }
                }
            }
        }

        double hAddV = 0.0;
        for(ui vzero = 0; vzero <= binNodes[0].size(); vzero++) {
            if(q<=h+vzero) break;
            ui vNzero = q-h-vzero;
            double chooseVNzeroNodes = 0;
            for(ui miss = 1; miss <= s-missEdges; miss++) {
                chooseVNzeroNodes += dp[vNzero][miss];
            }
            hAddV += chooseVNzeroNodes * CN[binNodes[0].size()][vzero];
        }
        for(ui u : H) {
            answers[u] += hAddV;
        }
      
        //addP
        for(ui xx = 0; xx <= s-missEdges; xx++) {
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.8\n");fflush(stdout);
}
#endif  
            if(binNodes[xx].size() == 0) continue;
            for(ui i = 0; i <= sumNodes; i++) {
                for(ui k = 0; k <= s-missEdges; k++) {
                    dp[i][k] = 0;
                }
            }
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.3\n");fflush(stdout);
}
#endif  
            for(ui miss = 1, i = 0; miss <= s-missEdges; miss++) {
                for(ui x = 0; x < binNodes[miss].size(); x++) {
                    i++;
                    for(ui j = i; j >= 1; j--) {
                        for(ui k = miss; k <= s-missEdges; k++) {
                            dp[j][k] += dp[j-1][ k-miss ];
                        }
                    }
                }
            }
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.4\n");fflush(stdout);
}
#endif  
            double hAddV = 0.0;
            for(ui vzero = 0; vzero <= binNodes[0].size() - (xx==0?1:0); vzero++) {
                if(q <= 1+h+vzero) break;
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.5\n");fflush(stdout);
}
#endif  
                ui vNzero = q-1-h-vzero;
                double chooseVNzeroNodes = 0;
                for(ui miss = 1; miss <= s-missEdges; miss++) {
                    chooseVNzeroNodes += dp[vNzero][miss];
                }
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.6 %u %.0f %u %u\n", vNzero,  chooseVNzeroNodes, binNodes[0].size()- (xx==0?1:0), vzero);fflush(stdout);
}
#endif  
                hAddV += chooseVNzeroNodes * CN[binNodes[0].size()- (xx==0?1:0)][vzero];
            }
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.7, haddV %.0f\n", hAddV);fflush(stdout);
}
#endif  
            for(ui u : binNodes[xx]) answers[u] += hAddV;
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.77 %u\n", xx);fflush(stdout);
}
#endif  
        }
    };

    if(sz == 0) {
        updateAns();
        
        return;
    }
    if(sz == 1) {
// P.push_back(C[0]);
binNodes[h-neiInH[C[0]]].push_back(C[0]);
    updateAns();
// P.pop_back();
binNodes[h-neiInH[C[0]]].pop_back();

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
#ifdef DEBUG2
if(uu==uuu) {
    printf("h1.5\n");fflush(stdout);
}
#endif
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
printf("pivot %u\n", pivot);
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
#ifdef DEBUG2
if(uu==uuu) {
    printf("h2\n");fflush(stdout);
}
#endif
    // if(neiInH[pivot] == h) {
    //     listing(deep+1, nxtC, newSz, p+1, h, missEdges);
    // }
    // else {
// P.push_back(pivot);
binNodes[h-neiInH[pivot]].push_back(pivot);
    listing(deep+1, nxtC, newSz, p+1, h, missEdges);
// P.pop_back();
binNodes[h-neiInH[pivot]].pop_back();
    // }

    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui v = sg.pEdge[i]; 
        // neiInH[v]--;
        level[v] = deep;
    }

    // level[pivot] = deep-1;
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

H.push_back(u);
        listing(deep+1, nxtC, newSz, p, h+1, newMissEdges);
H.pop_back();

        for(ui j = 0; j < newSz; j++) level[nxtC[j]] = deep;
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui v = sg.pEdge[j]; 
            if(level[v] != deep) continue;
            neiInH[v]--;
        }
    }

    freeMem(candSise);

}

sdcLocal::sdcLocal(Graph && g, ui s, ui q):g(g), s(s), q(q) {
    ui maxDepth = g.coreNumber * g.coreNumber;
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

    answers.resize(g.n);

    H.resize(maxDepth);
    H.clear();
    dp.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        dp[i].resize(s + 1);
    }
    dp[0][0] = 1;

    maxD = maxDepth;
    maxD2 = q + 1;
    computeC();

    binNodes.resize(s + 1);

    printf("sdcLocal::sdcVLocal\n");
}