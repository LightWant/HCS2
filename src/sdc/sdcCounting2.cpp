#include "sdcCounting2.h"
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
    printf("sdcCounting2.cpp::run");

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
#ifdef DDEBUG
for(ui i = 0; i < sz; i++) if(neiInPH[C[i]] != 0) printf("no 0: %u\n", C[i]);
// for(auto v:C) assert(neiInP[v] == 0);
#endif 
        for(ui i = 0; i < edC1; i++) neiInPH[C[i]] = 1;
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
double preAns = answers[q];
#endif

        listing(1, C, edC1, sz, 0, 1, 0);

        for(ui i = 0; i < sz; i++) neiInPH[C[i]] = 0;
        for(ui i = 0; i < sz; i++) level[C[i]] = 0;
        level[u] = 0;

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
std::vector<ui> cnti(g.n);
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        cnt++;
        // print(i);
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

void sdcCounting::pivoter(ui deep, ui * C, ui sz, ui p, ui h) {
#ifdef DDEBUG
printf("pdeep %u, p %u, h %u\nC:", deep, p, h);
for(ui i = 0; i < sz; i++) printf("%u ", C[i]); printf("\n");
#endif

    if(h == Q-1) {
        answers[Q-1] += 1;
        answers[Q] += p + sz;
#ifdef DDEBUG
printf("ans += %u, pivoter1\n", p + sz);
#endif
        return;
    }
    if(h+p+sz < q) return;

    auto updateAns = [&]() {
        // if(q>=h)
            // answers[q] += CN[p][q-h];
        for(ui i = std::max(h, q); i <= Q && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
#ifdef DDEBUG
if(i==Q)
printf("ans += C[%u][%u] %.0f, pivoter2\n", p, i-h, CN[p][i-h]);
#endif
        }
    };

    if(sz == 0) {
        updateAns(); return;
    }
    if(sz == 1) {
        p += 1; updateAns(); return;
    }
    if(sz == 2) {
        if(sg.deg[deep][C[0]] > sg.pIdx[C[0]]) {
            p += 2; updateAns();
        }
        else {
            p += 1; updateAns();
            p -= 1; h += 1; updateAns();
        }
        return;
    }

    ui pivot = C[0], pivotDeg = 0, pivotIndex = 0, num = 0;
    for(ui i = 0; i < sz; i++) {
        ui tmp = sg.deg[deep][C[i]] - sg.pIdx[C[i]];

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }
#ifdef DDEBUG
printf("pivot %u, pivotdeg %u, num %u\n", pivot, pivotDeg, num);
#endif
    if(pivotDeg+1 == sz && num == sz) {
        p += sz; updateAns(); return;
    }

    ui * nxtC = C + sz;
    ui newSz = 0;

    auto updateSG = [&]() {
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

    newSz = 0;
    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui w = sg.pEdge[i];
        nxtC[newSz++] = w;
        level[w] = deep + 1;
    }

    updateSG();

    ui candSize = sz - pivotDeg - 1;
    ui * cand = allocMem(candSize);
    for(ui i = 0, j = 0; i < sz; i++) {
        ui v = C[i];
        if(level[v] != deep + 1 && v != pivot) {
            cand[j++] = v;
        }
    }

    pivoter(deep + 1, nxtC, newSz, p + 1, h);
    for(ui i = 0; i < newSz; i++) level[nxtC[i]] = deep;

    for(ui i = 0; i < candSize; i++) {
        ui u = cand[i];
        level[u] = deep - 1;

        newSz = 0;
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui w = sg.pEdge[j];
            if(level[w] == deep - 1) continue;
            level[w] = deep + 1;
            nxtC[newSz++] = w;
        }
        updateSG();

        pivoter(deep + 1, nxtC, newSz, p, h + 1);

        for(ui i = 0; i < newSz; i++) level[nxtC[i]] = deep;
    }

    freeMem(candSize);
}

//C, cm, sz, p, h, missEdges, neiInPH, level
void sdcCounting::listing(ui deep, ui * C, ui cm, ui sz, ui p, ui h, ui missEdges) {
#ifdef DDEBUG
printf("       deep %u, p %u, h %u, missE %u\n", deep, p, h, missEdges);
printf("C_cm:");
for(ui i = 0; i < cm; i++) printf("%u ", C[i]); 
printf("C_sz:");
for(ui i = cm; i < sz; i++) printf("%u ", C[i]); 
printf("neiInP:");
for(ui i = 0; i < sz; i++) printf("%u ", neiInPH[C[i]]); printf("\n");
// printf("P:"); for(auto v:P) printf("%u ", v); printf("\n");

printf("\n");
#endif

    if(h == Q-1) {
// #ifdef DDEBUG
// printf("ans+%u\n", sz+p);
// #endif
        answers[Q-1] += 1;
        answers[Q] += sz + p;
        return;
    }
    if(h+p+sz < q) return;

    if(sz == 0) {
#ifdef DDEBUG
printf("ans_1+%.0f\n", CN[p][q-h]);
#endif
        // if(p >= q-h)
        //     answers[q] += CN[p][q-h];
        for(ui i = h>=q ? 0 : q-h; i+h <= Q && i <= p; i++) {
            answers[i + h] += CN[p][i];
        }
        return;
    }
    if(sz == 1) {
        p += 1;

        // if(p >= q-h)
        //     answers[q] += CN[p][q-h];
        for(ui i = h>=q ? 0 : q-h; i+h <= Q && i <= p; i++) {
            answers[i + h] += CN[p][i];
#ifdef DDEBUG
if(i+h==Q)
printf("ans_2+%.0f, C[%u][%u]\n", CN[p][Q-h], p, Q-h);
#endif
        }
        return;
    }

    if(missEdges == s) {
#ifdef DDEBUG
printf("pivoter");
#endif
        for(ui i = cm; i < sz; i++) level[C[i]] = deep - 1;
        for(ui j = 0; j < cm; j++) {
            ui u = C[j];
            ui & ed = sg.deg[deep][u];
            for(ui l = sg.pIdx[u]; l < ed; ) {
                ui w = sg.pEdge[l];
                if(level[w] != deep) {
                    std::swap(sg.pEdge[l], sg.pEdge[--ed]);
                }
                else l++;
            }
        }
        pivoter(deep, C, cm, p, h);
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
    auto addPNode = [&](ui u, ui & newCm, ui & newSz) {
        for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
            ui v = sg.pEdge[i]; 
            neiInPH[v]++;
            level[v] = deep+1;
        }

        for(ui i = 0; i < cm; i++) {
            ui v = C[i];
            if(level[v] == deep+1) {
                nxtC[newCm++] = v;
            }
        }
        
        newSz = newCm;
        for(ui i = cm; i < sz; i++) {
            ui v = C[i];
            if(level[v] == deep+1) {
                nxtC[newSz++] = v;
            }
        }

        updateSG(newSz);
    };

    auto addHNode = [&](ui u, ui & newCm, ui & newSz, ui newMissEdges) {
        for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
            ui v = sg.pEdge[i]; 
            if(level[v] != deep) continue;
            neiInPH[v]++;
            level[v] = deep+2;
        }

        for(ui i = 0; i < cm; i++) {
            ui v = C[i];
            if(level[v] != deep+2) continue;

            if(newMissEdges + p+h+1-neiInPH[v] <= s) {
                level[v] = deep+1;
                nxtC[newCm++] = v;
            }
            else level[v] = deep;
        }

        newSz = newCm;

        for(ui i = 0; i < cm; i++) {
            ui v = C[i];
            if(level[v] != deep) continue;

            if(newMissEdges + p+h+1-neiInPH[v] <= s) {
                level[v] = deep+1;
                nxtC[newSz++] = v;
            }
        }

        for(ui i = cm; i < sz; i++) {
            ui v = C[i];
            if(level[v] < deep) continue;
// printf("szzz %u, %u %u %u\n", v, newMissEdges, p+h+1, neiInPH[v]);
            if(newMissEdges + p+h+1-neiInPH[v] <= s) {
                level[v] = deep+1;
                nxtC[newSz++] = v;
            }
            else if(level[v] == deep+2) level[v] = deep;
        }
// for(ui i = 0; i < sz; i++) {
//     bool in = false;
//     for(ui j = 0; j < newSz; j++) {
//         if(nxtC[j] == C[i]) {
//             in = true; break;
//         }
//     }

//     if(in) assert(level[C[i]] == deep + 1);
//     else assert(level[C[i]] <= deep); 
// }
        updateSG(newSz);
    };
    auto rmHNode = [&](ui u, ui newSz) {
        for(ui i = 0; i < newSz; i++) level[nxtC[i]] = deep;
// for(ui i = 0; i < sz; i++) assert(level[C[i]] <= deep);
        for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
            ui v = sg.pEdge[i]; 
            if(level[v] != deep) continue;
            neiInPH[v]--;
        }
    };

    if(cm > 0) {
        ui pivot = C[0], pivotDeg = 0;
        for(ui i = 0; i < cm; i++) {
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
        ui newCm = 0, newSz = 0;
        addPNode(pivot, newCm, newSz);

        ui candSise = sz - pivotDeg - 1;
        ui * cand = allocMem(candSise);
        for(ui i = 0, j = 0; i < sz; i++) {
            ui v = C[i];
            if(level[v] != deep + 1 && v != pivot) {
                cand[j++] = v;
// assert(j <= candSise);
            }
        }

        listing(deep+1, nxtC, newCm, newSz, p+1, h, missEdges);

        for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
            ui v = sg.pEdge[i]; 
            neiInPH[v]--;
            level[v] = deep;
        }
#ifdef DDEBUG
printf("       After deep %u, p %u, h %u, missE %u\n", deep, p, h, missEdges);
printf("C_cm:");
for(ui i = 0; i < cm; i++) printf("%u ", C[i]); 
printf("C_sz:");
for(ui i = cm; i < sz; i++) printf("%u ", C[i]); 
printf("neiInP:");
for(ui i = 0; i < sz; i++) printf("%u ", neiInPH[C[i]]); printf("\n");
printf("cands:");
for(ui i = 0; i < candSise; i++) printf("%u ", cand[i]); printf("\n");

printf("\n");
#endif

        for(ui i = 0; i < candSise; i++) {
            ui u = cand[i];
#ifdef DDEBUG
printf("  deep %u, cand %u, neiIPH %u\n", deep, u, neiInPH[u]);
printf("cands:");
for(ui i = 0; i < candSise; i++) printf("%u ", cand[i]); printf("\n");

#endif
            level[u] = deep-1;
            ui newCm = 0, newSz = 0, newMissEdges = missEdges + p+h-neiInPH[u];
            addHNode(u, newCm, newSz, newMissEdges);

            listing(deep+1, nxtC, newCm, newSz, p, h+1, newMissEdges);

            rmHNode(u, newSz);
        }

        freeMem(candSise);
    }
    else {

        for(ui i = h>=q ? 0 : q-h; i+h <= Q && i <= p; i++) {
            answers[i + h] += CN[p][i];
#ifdef DDEBUG
if(i+h==Q)
printf("ans_3+%.0f\n", CN[p][Q-h]);
#endif
        }
        for(ui i = 0; i < sz; i++) {
            ui u = C[i];
#ifdef DDEBUG
printf("  edeep %u, cand %u, neiIPH %u\n", deep, u, neiInPH[u]);
#endif
            level[u] = deep-1;
            ui newCm = 0, newSz = 0, newMissEdges = missEdges + p+h-neiInPH[u];
            addHNode(u, newCm, newSz, newMissEdges);

            listing(deep+1, nxtC, newCm, newSz, p, h+1, newMissEdges);

            rmHNode(u, newSz);
        }
    }
}

sdcCounting::sdcCounting(Graph && g, ui s, ui q, ui Q):g(g), s(s), q(q), Q(Q) {
    ui maxDepth = g.coreNumber + s + 1;
    ui maxCSize = std::min(g.coreNumber * g.coreNumber, g.n);
    candBuffer = new ui[maxDepth * maxCSize];

    initBuffer(maxDepth * maxCSize);

    sg.pIdx.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(maxDepth);
    for(ui i = 0; i < maxDepth; i++) {
        sg.deg[i].resize(g.n);
    }

    neiInPH.resize(g.n);

    level.resize(g.n);

    answers.resize(Q + 1);

    maxD = maxDepth;
    maxD2 = Q + 1;
    computeC();

    printf("sdcCounting2::sdcCounting\n");
}

sdcCounting::~sdcCounting() {
    if(bf3 != nullptr) delete [] bf3;
    if(CN != nullptr) delete [] CN;
    if(candBuffer != nullptr) delete [] candBuffer;
    if(memBuffer != nullptr) delete [] memBuffer;
}

