#include "maximalSDC.h"

#include <queue>
#include <utility>

#define BASELINE
#define DDEBUG

#ifdef DDEBUG
#include <iostream>
#endif

//s <= q-2. 1.dinamiter2 2. 连通
//common nei q-2-(s), thus at least q-s-2 common neighbors
ull msdc::run() {
#ifdef DDEBUG
g.print();
#endif
    printf("maximalSDC.cpp\n");

    g.initHash();
    printf("init Hash\n");

    // std::vector<ui> toDelete;
    // toDelete.clear();

    using Pair = std::pair<ui, ui>;
    std::queue<ui> que;

    for(ui u = 0; u < g.n; u++) {
        ui edC = 0;
        ui stX = 0, edX = 0;
#ifdef DDEBUG
std::cout<<"    start "<<u<<' '<<answer<<std::endl;
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
            if(deg[v] >= q-s-2) C.changeTo(v, edC++);
            deg[v] = 0;
        }

#ifdef DDDEBUG
printf("C 1-hop:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
#endif
        //X in 1-hop nei
        for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
            ui v = g.pEdge[i];
            ui d = 0;
            for(ui i = 0, ed = edC; i < ed; i++) {
                if(g.connectHash(u, C[i])) d++;
            }
            if(d >= q-s-2) X.changeTo(v, edX++);
        }
        
#ifdef DDDEBUG
printf("X 1-hop:"); for(ui i = 0; i < edX; i++) printf("%u ", X[i]); printf("\n");
#endif
        //2-hop
        ui edC1 = edC, edX1 = edX;
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connectHash(u, w) || w == u) continue;

                if((++deg[w]) == q-s-2 || (q==s+2 && deg[w]==1)) {
                    if(w < u) X.changeTo(w, edX++);
                    else C.changeTo(w, edC++);
                }
            }
        }
        for(ui i = 0; i < edC1; i++) {
            ui v = C[i];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                deg[w] = 0;
            }
        }
#ifdef DDDEBUG
printf("C 2-hop:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("X 2-hop:"); for(ui i = 0; i < edX; i++) printf("%u ", X[i]); printf("\n");
#endif
        //build sub-graph g
        auto buildV = [&](ui v) {
            sg.pIdx[v] = sg.pIdx2[v] = g.pIdx[v];

            if(g.connectHash(v, u)) sg.pEdge[sg.pIdx2[v]++] = u;
            for(ui i = 0; i < edC; i++) {
                if(g.connectHash(v, C[i])) sg.pEdge[sg.pIdx2[v]++] = C[i];
            }
            for(ui i = 0; i < edX; i++) {
                if(g.connectHash(v, X[i])) sg.pEdge[sg.pIdx2[v]++] = X[i];
            }
        };
        for(ui i = 0; i < edC; i++) buildV(C[i]);
        for(ui i = stX; i < edX; i++) buildV(X[i]);
        buildV(u);

        for(ui i = 0; i < edC1; i++) neiInP[C[i]] = 1;
        for(ui i = 0; i < edX1; i++) neiInP[X[i]] = 1;
        P.push_back(u);
        bkPivot(1, stX, edX, edC, 0);
        P.pop_back();
        for(ui i = 0; i < edC; i++) neiInP[C[i]] = 0;
        for(ui i = 0; i < edX; i++) neiInP[X[i]] = 0;
    }

#ifdef BASELINE
#define BASELINES 2
#define BASELINEQ 5
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto getMissEdges = [&](ui x) {
    ui TotalD = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
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
    if(sz < BASELINEQ) return false;
    return getMissEdges(x) <= BASELINES;
};

ui cnt = 0;
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        //isMaximal
        bool isMaximal = true;
        for(ui u = 0; u < g.n; u++) {
            if((1<<u) & i) continue;
            if(check(i | (1<<u))) {
                isMaximal = false;
                break;
            }
        }
        if(isMaximal) {
            cnt++;
            print(i);
        }
    }
}
printf("cnt:%u\n", cnt);
#endif

    return answer;
}

void msdc::bkPivot(ui deep, ui stX, ui edX, ui edC, ui missedEdges) {
#ifdef DDEBUG
printf("       deep %u\n", deep);
printf("C:");
for(ui i = 0; i < edC; i++) printf("%u ", C[i]); 
printf("neiInPC:");
for(ui i = 0; i < edC; i++) printf("%u ", neiInP[C[i]]); printf("\n");
printf("X:");
for(ui i = stX; i < edX; i++) printf("%u ", X[i]); 
printf("neiInPX:");
for(ui i = stX; i < edX; i++) printf("%u ", neiInP[X[i]]); printf("\n");
printf("P:"); for(auto v:P) printf("%u ", v); printf("\n");
printf("missE %u\n", missedEdges);
printf("\n");
#endif
    if(P.size() + edC < q) return;
    if(edC == 0) {
        if(stX == edX) {
#ifdef DDEBUG
printf("ans normal\n");
#endif
            answer++;
        }
#ifdef DDEBUG
printf("X not empty\n");
#endif
        return;
    }

    //init deg from the pre level
    if(sg.deg[deep].size() == 0) sg.deg[deep].resize(g.n);
    auto computeDegInC = [&](ui v) {
        ui costC = edC*2;
        ui costAdj = sg.degree2(v);
        ui d = 0;
        if(costC < costAdj) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) d++;
        }
        else  {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                ui w = sg.pEdge[i];
                if(C.isIn(w, 0, edC)) d++;
            }
        }
        sg.deg[deep][v] = d;
    };
    for(ui i = 0; i < edC; i++) computeDegInC(C[i]);
    for(ui i = 0; i < P.size(); i++) computeDegInC(P[i]);
    for(ui i = stX; i < edX; i++) computeDegInC(X[i]);

    auto hasUniversalX = [&]() {
        for(ui i = stX; i < edX; i++) {
            if(neiInP[X[i]] < P.size()) continue;
            if(sg.deg[deep][X[i]] == edC) return true;
        }
        return false;
    };
    if(hasUniversalX()) {
#ifdef DDDEBUG
printf("Universial\n");
#endif
        return;
    }

    //is P+C a k-plex
    ui totalMissedEdges = missedEdges;
    for(ui i = 0; i < edC; i++) {
        ui v = C[i];
        totalMissedEdges += P.size() - neiInP[v];
    }
    if(totalMissedEdges <= s) {
        ui missEdgesInC = 0;
        for(ui i = 0; i < edC; i++) {
            missEdgesInC += edC-sg.deg[deep][C[i]]-1;
        }
        totalMissedEdges += missEdgesInC / 2;
    }

    if(totalMissedEdges <= s) {//P+C is a sdc, check X
#ifdef DDEBUG
printf("totalMissedEdges %u\n", totalMissedEdges);
#endif
        for(ui i = stX; i < edX; i++) {
            ui v = X[i];

            ui medges = totalMissedEdges + 
                        P.size() - neiInP[v] +
                        edC - sg.deg[deep][v];
            if(medges > s) continue;
            return;
        }
#ifdef DDEBUG
printf("ans P+C\n");
#endif
        // bkPivot(deep + 1, edX, edX, 0, totalMissedEdges);
        answer++;
        return;
    }

    //find pivot in C and X with max degree in C
    ui maxD = 0, pv = C[0];
    bool inC = true;

    for(ui i = 0; i < edC; i++) {
        ui v = C[i];
        ui d = sg.deg[deep][v];

        if(d > maxD) {
            maxD = d;
            pv = v;
        }
    }
    for(ui i = stX; i < edX; i++) {
        ui v = X[i];
        ui d = sg.deg[deep][v];
        
        if(d > maxD) {
            maxD = d;
            pv = v;
            inC = false;
        }
    }
#ifdef DDEBUG
printf("pv %u\n", pv);
#endif
    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;

    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        // if(g.connectHash(pv, v)) {
        //     if(neiInP[v] != P.size()) cand[candSize++] = v;
        // }
        // else cand[candSize++] = v;
        cand[candSize++] = v;
    }
#ifdef DDEBUG
printf("cands:");
for(ui i = 0; i < candSize; i++) printf("%u ", cand[i]);printf("\n");
#endif
    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];
#ifdef DDEBUG
printf("        deep %u, i %u, v %u\n", deep, i, v);
printf("C:");
for(ui i = 0; i < edC; i++) printf("%u ", C[i]); 
printf("neiInPC:");
for(ui i = 0; i < edC; i++) printf("%u ", neiInP[C[i]]); printf("\n");
printf("X:");
for(ui i = stX; i < edX; i++) printf("%u ", X[i]);
printf("neiInPX:");
for(ui i = stX; i < edX; i++) printf("%u ", neiInP[X[i]]); printf("\n");
printf("P:"); for(auto v:P) printf("%u ", v); printf("\n");
printf("missE %u\n", missedEdges);
#endif
        ui newMissedEdges = missedEdges + P.size() - neiInP[v];

#ifdef DDEBUG
printf("newMissedEdges %u\n", newMissedEdges);
#endif
        if(newMissedEdges > s) {
            C.changeTo(v, --edC);
            X.changeTo(v, edX++);
            continue;
        }

        P.push_back(v);
        C.changeTo(v, --edC);

        //build tmpC, tmpStX and tmpEdX for the next level
        ui tmpC = 0;
        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(v, C[j])) {
// #ifdef DDEBUG
// printf("Nei %u, newM %u, P %u neiInP[C[j]] %u %u\n", C[j], newMissedEdges, 
// P.size(), neiInP[C[j]],
// newMissedEdges + P.size() - neiInP[C[j]]);
// #endif
                if(newMissedEdges + P.size() - neiInP[C[j]] - 1 <= s) {
                    neiInP[C[j]]++;
                    C.changeToByPos(j, tmpC++);
                }
            }
            else {
// #ifdef DDEBUG
// printf("nonNei %u, %u, %u %u %u\n", C[j], newMissedEdges, P.size(), neiInP[C[j]],
// newMissedEdges + P.size() - neiInP[C[j]] + 1);
// #endif
                if(newMissedEdges + P.size() - neiInP[C[j]] <= s) {
                    C.changeToByPos(j, tmpC++);
                }
            }
        }
        ui tmpStX = edX;
        for(ui j = stX; j < tmpStX; ) {
            if(g.connectHash(v, X[j])) {
                if(newMissedEdges + P.size() - neiInP[X[j]] - 1 <= s) {
                    neiInP[X[j]]++; 
                    X.changeToByPos(j, --tmpStX);
                }
                else j++;
            }
            else {
                if(newMissedEdges + P.size() - neiInP[X[j]] <= s) {
                    X.changeToByPos(j, --tmpStX);
                }
                else j++;
            }
        }

        bkPivot(deep + 1, tmpStX, edX, tmpC, newMissedEdges);

        for(ui j = tmpStX; j < edX; j++) {
            if(g.connectHash(v, X[j])) neiInP[X[j]]--;
        }
        for(ui j = 0; j < tmpC; j++) {
            if(g.connectHash(v, C[j])) neiInP[C[j]]--;
        }
        P.pop_back();
        X.changeTo(v, edX++);
    }

    if(candSize > 0) {
        for(ui i = candSize - 1; i > 0; i--) {
            ui v = cand[i];
            X.changeTo(v, --edX);
        }
        X.changeTo(cand[0], --edX);
    }
}

// forEachCommonElementDo(g.pEdge + g.pIdx2[u], g.degree2(u),
//                 g.pEdge + g.pIdx[v], g.degree(v),
//                 [&](ui w) {
//                     if(deg[w] >= q-s-1 && (--deg[w]) < q-s-1) {
//                         que.push({w, deg[w]});
//                     }
//                 };
//             );