#include "sdcLocal.h"
#include <queue>
#include <algorithm>
#include <utility>

// #define BASELINE
// #define DDEBUG

#ifdef DDEBUG
#include <iostream>
#endif

void sdcLocal::run() {
    printf("sdcELocal.cpp::run\n");

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

    auto getEdge = [&](ui a, ui b) {
        if(a > b) std::swap(a, b);

        auto st = pOEdge.begin() + pOIdx[a];
        auto ed = pOEdge.begin() + pOIdx[a + 1];
        ui idx = std::lower_bound(st, ed, b) - pOEdge.begin(); 

        return idx;
    };

    if(h == q-1) {
        for(ui i = 0; i < sz; i++) {
            ui u = C[i];
            for(auto v : H) if(g.connectHash(u, v)) {
                ui idx = getEdge(u, v);
                answers[idx] += 1;
#ifdef DDEBUG
printf("edge %u-%u, add 1 (h-C)\n", u, v, 1);
#endif
                
            }
        }
        ui cntOfNodesSmall = 0;
        for(ui i = 0; i <= s-missEdges; i++) {
            cntOfNodesSmall += binNodes[i].size();

            for(auto u : binNodes[i]) {
                bool ok = false;
                for(auto v : H) if(g.connectHash(u, v)) {
                    ui idx = getEdge(u, v);
                    answers[idx] += 1;
#ifdef DDEBUG
printf("edge %u-%u, add 1 (ph)\n", u, v, 1);
#endif
                }
            }
        }
        for(ui i = 0; i < H.size(); i++) 

        if(cntOfNodesSmall+sz > 0)
        for(ui j = i + 1; j < H.size(); j++) if(g.connectHash(H[i], H[j])) {
            ui idx = getEdge(H[i], H[j]);
            answers[idx] += sz + cntOfNodesSmall;
#ifdef DDEBUG
printf("edge %u-%u, add %u (hh)\n", H[i], H[j], sz + cntOfNodesSmall);
#endif
        }
        return;
    }
    if(h+p+sz < q) return;

    auto updateAns = [&]() {
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
        // for(each kind of edge) {
        //     get count
        //     for each edge in this kind:
        //         add += count
        // }
        ui sumNodes = 0;
        for(ui i = 0; i <= s-missEdges; i++) sumNodes += binNodes[i].size();
        
        //pp
        for(ui i = 0; i <= s-missEdges; i++) {
            if(binNodes[i].size() == 0) continue;
            for(ui j = i; j <= s-missEdges; j++) if(i+j <= s-missEdges) {
                if(binNodes[j].size() == 0) continue;
                if(j == i && binNodes[i].size() == 1) continue;
                //dp[x][w] = dp[x][w] + dp[x-1][w-sid]
                
                ui edgesToChoose = s-missEdges-i-j;
                ui nodesToChoose = q-h-2;

                for(ui j = 1; j <= nodesToChoose; j++) {
                    for(ui k = 0; k <= edgesToChoose; k++) {
                        dp[j][k] = 0;
                    }
                }
                for(ui t = 0; t <= edgesToChoose; t++) {
                    ui itrs = binNodes[t].size();
                    if(t == i) itrs--;
                    if(t == j) itrs--;

                    for(ui x = 1; x <= itrs; x++) {
                        //choose x nodes from itrs nodes
                        for(ui y = nodesToChoose; y >= 1; y--) {
                            for(ui z = t; z <= edgesToChoose; z++) {
                                dp[y][z] += dp[y - 1][z - t];
                            }
                        }
                    }
                }

                ui cntOfWays = 0;
                for(ui x = 0; x <= edgesToChoose; x++) {
                    cntOfWays += dp[nodesToChoose][x];
                }
#ifdef DDEBUG
printf("p %u, p %u, cnt %u\n", i, j, cntOfWays);
#endif
                if(i != j)
                for(auto u : binNodes[i]) {
                    for(auto v : binNodes[j]) {
                        ui idx = getEdge(u, v);
                        answers[idx] += cntOfWays;
#ifdef DDEBUG
printf("edge %u-%u, add %u\n", u, v, cntOfWays);
#endif
                    }
                }
                else {
                    for(ui x = 0; x < binNodes[i].size(); x++) {
                        for(ui y = x + 1; y < binNodes[i].size(); y++) {
                            ui idx = getEdge(binNodes[i][x], binNodes[i][y]);
                            answers[idx] += cntOfWays;
#ifdef DDEBUG
printf("edge %u-%u, add %u\n", binNodes[i][x], binNodes[i][y], cntOfWays);
#endif
                        }
                    }
                }
            }
        }

        //ph
        for(ui i = 0; i <= s-missEdges; i++) {
            if(binNodes[i].size() == 0) continue;

            ui edgesToChoose = s-missEdges-i;
            ui nodesToChoose = q-h-1;

            for(ui j = 1; j <= nodesToChoose; j++) {
                for(ui k = 0; k <= edgesToChoose; k++) {
                    dp[j][k] = 0;
                }
            }    

            for(ui t = 0; t <= edgesToChoose; t++) {
                ui itrs = binNodes[t].size();
                if(t == i) itrs--;

                for(ui x = 1; x <= itrs; x++) {
                    //choose x nodes from itrs nodes
                    for(ui y = nodesToChoose; y >= 1; y--) {
                        for(ui z = t; z <= edgesToChoose; z++) {
                            dp[y][z] += dp[y - 1][z - t];
                        }
                    }
                }
            }

            ui cntOfWays = 0;
            for(ui x = 0; x <= edgesToChoose; x++) {
                cntOfWays += dp[nodesToChoose][x];
            }
#ifdef DDEBUG
printf("p %u, h, %u\n", i, cntOfWays);
#endif
            for(auto u : binNodes[i]) {
                for(auto v : H) if(g.connectHash(u, v)) {
                    ui idx = getEdge(u, v);
                    answers[idx] += cntOfWays;
#ifdef DDEBUG
printf("edge %u-%u, add %u\n", u, v, cntOfWays);
#endif
                }
            }
        }

        //hh
        ui edgesToChoose = s-missEdges;
        ui nodesToChoose = q-h;

        for(ui j = 1; j <= nodesToChoose; j++) {
            for(ui k = 0; k <= edgesToChoose; k++) {
                dp[j][k] = 0;
            }
        }

        for(ui t = 0; t <= edgesToChoose; t++) {
            ui itrs = binNodes[t].size();

            for(ui x = 1; x <= itrs; x++) {
                //choose x nodes from itrs nodes
                for(ui y = nodesToChoose; y >= 1; y--) {
                    for(ui z = t; z <= edgesToChoose; z++) {
                        dp[y][z] += dp[y - 1][z - t];
// #ifdef DDEBUG
// printf("t %u, dp[y %u][z %u] += dp[y-x %u][z-t*x %u] %.0f\n", 
//     t, y, z, y-x, z-t*x, dp[y - x][z - t*x]);
// #endif
                    }
                }
            }
        }

        ui cntOfWays = 0;
        for(ui x = 0; x <= edgesToChoose; x++) {
            cntOfWays += dp[nodesToChoose][x];
        }
#ifdef DDEBUG
printf("hh %u\n", cntOfWays);
// for(ui x = 0; x <= edgesToChoose; x++) {
//     printf("s%u-%.0f ", x, dp[nodesToChoose][x]);
// }
// printf("\n");
#endif
        for(ui i = 0; i < H.size(); i++) {
            ui u = H[i];
            for(ui j = i + 1; j < H.size(); j++) {
                ui v = H[j];
                if(g.connectHash(u, v)) {
                    ui idx = getEdge(u, v);
#ifdef DDEBUG
printf("edge %u-%u, add %u\n", u, v, cntOfWays);
#endif
                    answers[idx] += cntOfWays;
                }
            }
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

    answers.resize(g.m);

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

    pOEdge.resize(g.m / 2);
    pOIdx.resize(g.n + 1);
    for(ui u = 0; u < g.n; u++) {
        pOIdx[u + 1] = pOIdx[u];
        for(ui j = g.pIdx2[u]; j < g.pIdx[u + 1]; j++) {
            pOEdge[pOIdx[u+1]++] = g.pEdge[j];
        }
    }

    //only when g.edges is sorted
    reEegeId.resize(g.m/2);
    for(ui i = 0; i < g.m/2; i++) {
        ui u = g.edges[i].first;
        ui v = g.edges[i].second;
        u = g.mp2[u];
        v = g.mp2[v];
        if(u > v) std::swap(u, v);

        auto st = pOEdge.begin() + pOIdx[u];
        auto ed = pOEdge.begin() + pOIdx[u + 1];
        ui idx = std::lower_bound(st, ed, v) - pOEdge.begin(); 

        reEegeId[i] = idx;
    }

    printf("sdcLocal::sdcELocal\n");fflush(stdout);
}