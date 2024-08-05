#include "d2k.h"
// #define DDDEBUG
// #define BASELINE

ui d2k::run() {
#ifdef DDDEBUG
g.print();
#endif
    printf("d2ksginC\n");
    g.initHash();
    printf("init Hash\n"); fflush(stdout);

    std::vector<ui> toDelete(g.maxD);
    toDelete.clear();

    // sg.pDeepIdx.resize(g.n);
    // sg.pDeepIdx[0].resize(g.n);
    sg.deg.resize(g.n);
    sg.deg[0].resize(g.n);

    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<' '<<answer<<std::endl;
#endif
        //P is empty
        //X is 2-hop neighbors < u
        //C is 2-hop neighbors > u
        if(g.pIdx[u + 1] - g.pIdx2[u] + k < q) continue;
        ui edC = 0;
        ui stX = 0, edX = 0;
        ui edP = 0;
        // nn.clear();

        P.changeTo(u, edP++);
        nn.addNonNei(u, u);

        for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) X.changeTo(g.pEdge[i], edX++);
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) C.changeTo(g.pEdge[i], edC++);
#ifdef DDDEBUG
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = 0; i < edX; i++) printf("%u ", X[i]); printf("\n");
printf("edX %u, edC %u\n", edX, edC);
#endif
        //delete C
        while(true) {
            for(ui i = 0; i < edC; i++) {
                ui v = C[i];
                ui d = 0;
                if(g.degree(v) <= edC * logBasedOn2[g.degree(v)]) {
                    for(ui i = g.pIdx[v]; i < g.pIdx[v + 1]; i++) {
                        if(C.isIn(g.pEdge[i], 0, edC)) {
                            d++;
                            if(d + 2*k >= q) break;
                        }
                    }
                }
                else {
                    for(ui j = 0; j < edC; j++) {
                        ui w = C[j];
                        if(g.connectHash(v, w)) {
                            d++;
                            if(d + 2*k >= q) break;
                        }
                    }
                }

                if(d + 2*k < q) {//
                    // C.remove(v);
                    toDelete.push_back(v);
                }
            }

            if(toDelete.size() == 0) break;
            for(auto v : toDelete) C.changeTo(v, --edC);
            toDelete.clear();
        }
#ifdef DDDEBUG
printf("After delete C\n");
#endif
        //delete X
        for(ui i = 0; i < edX; i++) {
            ui v = X[i];
            ui d = 0;

            if(g.degree(v) <= edC * logBasedOn2[g.degree(v)]) {
                for(ui i = g.pIdx[v]; i < g.pIdx[v + 1]; i++) {
                    if(C.isIn(g.pEdge[i], 0, edC)) {
                        d++;
                        if(d + 2*k >= q) break;
                    }
                }
            }
            else {
                for(ui j = 0; j < edC; j++) {
                    ui w = C[j]; 
                    if(g.connectHash(v, w)) {
                        d++;
                        if(d + 2*k >= q) break;
                    }
                }
            }
            
            if(d + 2*k < q) X.changeTo(v, --edX);
        }
#ifdef DDDEBUG
printf("After Delete X\n");
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = 0; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
// printf("neighbors in P:\n");
// for(ui i = 0; i < edC; i++) printf("%u ", neighborsInP[C[i]]); printf("\n");
// for(ui i = 0; i < edX; i++) printf("%u ", neighborsInP[X[i]]); printf("\n");
#endif
        //two hop neighbors
        if(k > 1)
        for(ui j = 0, ed = edC; j < ed; j++) {
            ui v = C[j];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connectHash(u, w)) continue;
                if(C.isIn(w, 0, edC) || X.isIn(w, 0, edX) || w == u) continue; 
                ui d = 0;
                
                for(ui j = 0; j < ed; j++) {
                    ui x = C[j]; 
                    if(g.connectHash(w, x)) {
                        d++;

                        if(d + 2 * k >= q + 2) {
                            if(w > u) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in C\n", w, v);
#endif
                                C.changeTo(w, edC++);
                                nn.addNonNei(w, u);
                            }
                            else if(w < u) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in X\n", w, v);
#endif
                                X.changeTo(w, edX++);
                                nn.addNonNei(w, u);
                            }
                            break;
                        }
                    }
                }
            }
        }

        
#ifdef DDDEBUG
printf("After Add 2-hop neighbors\n");
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = 0; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
#endif
        //build sub-graph g
        auto buildV = [&](ui v) {
            sg.pIdx[v] = sg.pIdx2[v] = g.pIdx[v];
            for(ui i = g.pIdx[v]; i < g.pIdx[v + 1]; i++) {

                ui w = g.pEdge[i];

                if(C.isIn(w, 0, edC)) {
                    sg.pEdge[sg.pIdx2[v]++] = w;
                }
            }
        };
        
        for(ui i = 0; i < edC; i++) buildV(C[i]);
        for(ui i = 0; i < edX; i++) buildV(X[i]);
        for(ui i = 0; i < edC; i++) sg.deg[0][C[i]] = sg.pIdx2[C[i]] - sg.pIdx[C[i]];
        for(ui i = 0; i < edX; i++) sg.deg[0][X[i]] = sg.pIdx2[X[i]] - sg.pIdx[X[i]];
        buildV(u);
        sg.deg[0][u] = sg.pIdx2[u] - sg.pIdx[u];
        for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
            ui w = g.pEdge[i];
            if(C.isIn(w, 0, edC) || X.isIn(w, 0, edX)) {
                sg.pEdge[sg.pIdx2[w]++] = u;
            }
        }
        
#ifdef DDDEBUG
printf("sg:\n");
for(ui i = 0; i < edC; i++) {
    ui v = C[i];
    printf("%u:", v);
    for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
        ui w = sg.pEdge[j];
        printf("%u ", w);
    }printf("\n");
}printf("\n");
for(ui i = 0; i < edX; i++) {
    ui v = X[i];
    printf("%u:", v);
    for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
        ui w = sg.pEdge[j];
        printf("%u ", w);
    }
}printf("\n");
#endif

        bkPivot(1, stX, edX, edC, edP);

        for(ui i = 0; i < edC; i++) nn.clear(C[i]);
        for(ui i = 0; i < edP; i++) nn.clear(P[i]);
        for(ui i = stX; i < edX; i++) nn.clear(X[i]);
    }

#ifdef BASELINE
#define BASELINEK 2
#define BASELINEQ 3
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz < BASELINEQ) return false;

    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
        ui d = 0;
        for(ui v = 0; v < g.n; v++) if((1<<v) & x){
            if(g.connectHash(u, v)) d++;
        }
        if(d + BASELINEK < sz) return false; 
    }

    return true;
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

void d2k::bkPivot(ui deep, ui stX, ui edX, ui edC, ui edP) {
#ifdef DDDEBUG
printf("deep %u\n", deep);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = stX; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
#endif

    auto isMaximal = [&]() {
        for(ui i = 0; i < edC; i++) if(nn.getCntNonNei(C[i]) != edP) return false;
        for(ui i = stX; i < edX; i++) if(nn.getCntNonNei(C[i]) != edP) return false;
        return true;
    };
    if(edC == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        if(stX == edX && edP >= q) {
#ifdef DDDEBUG
printf("    P is ans!\n");
#endif
            answer++;
        }
        return;
    }

    if(edP + edC < q) return;

    auto hasUniversalX = [&]() {
        for(ui i = stX; i < edX; i++) {
            if(nn.getCntNonNei(X[i])) continue;
            if(sg.deg[deep - 1][X[i]] == edC) return true;
// ui d = 0;
// for(ui j = 0; j < edC; j++) {
//     if(g.connectHash(X[i], C[j])) d++;
// }
// assert(sg.deg[deep - 1][X[i]] == d);
            // ui cntNN = nn.getCntNonNei(X[i]);
            // if(cntNN + 1 + edC - sg.deg[deep - 1][X[i]] <= k) {
            //     bool canAdd = true;
            //     for(ui j = 0; j < cntNN; j++) {
            //         ui w = nn.buffer[X[i]*k + j];
            //         if(nn.getCntNonNei(w) + edC - sg.deg[deep-1][w] >= k) {
            //             canAdd = false; break;
            //         }
            //     }
            //     if(!canAdd) continue;

            //     if(canAdd) return true;
            // }
        }
        return false;
    };
    

    //init deg from the pre level
    if(sg.deg[deep].size() == 0) sg.deg[deep].resize(g.n);
    for(ui i = 0; i < edC; i++) sg.deg[deep][C[i]]=sg.deg[deep-1][C[i]];
    for(ui i = 0; i < edP; i++) sg.deg[deep][P[i]]=sg.deg[deep-1][P[i]];
    for(ui i = stX; i < edX; i++) sg.deg[deep][X[i]]=sg.deg[deep-1][X[i]];
#ifdef DDDEBUG
printf("CdeginC:");
for(ui i = 0; i < edC; i++) printf("%u-%u ", C[i], sg.deg[deep][C[i]]);
printf("\nPdeginC:");
for(ui i = 0; i < edP; i++) printf("%u-%u ", P[i], sg.deg[deep][P[i]]);
printf("\nXdeginC:");
for(ui i = stX; i < edX; i++) printf("%u-%u ", X[i], sg.deg[deep][X[i]]);
printf("\n");
#endif

    if(hasUniversalX()) {
#ifdef DDDEBUG
printf("Universial\n");
#endif
        return;
    }
    //is P+C a k-plex
    ui minDPC = g.n;
    //find pivot in C and X with max degree in C
    ui maxD = 0, pv = C[0];
    bool inC = true;

    for(ui i = 0; i < edC; i++) {
        ui v = C[i];
        ui d = sg.degreeDeep(deep, v);
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
        if(d > maxD) {
            maxD = d;
            pv = v;
        }
    }
    for(ui i = stX; i < edX; i++) {
        ui v = X[i];
        ui d = sg.degreeDeep(deep, v);
        
        if(d > maxD) {
            maxD = d;
            pv = v;
            inC = false;
        }
    }

    //compute minDPC
    if(minDPC + k >= edP + edC)
    for(ui i = 0; i < edP; i++) {
        ui v = P[i];
        ui d = sg.degreeDeep(deep, v);
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
        if(minDPC + k < edP + edC) break;
    }

    if(minDPC + k >= edP + edC) {//P+C is a k-plex, check X
#ifdef DDDEBUG
printf(" P+C!; minDPC %u\n", minDPC);
#endif
        bool ok = true;

        for(ui i = stX; i < edX; i++) {
            ui u = X[i];
            bool canAdd = true;

            //check if P+C+u is a plex
            ui nnU = nn.getCntNonNei(u) + 1;

            if(nnU <= k)
            for(ui j = 0; j < edC; j++) {
                ui v = C[j];

                if(!g.connectHash(u, v)) {
                    nnU++;
                    if(nnU > k) {
                        canAdd = false;
                        break;
                    }
                    if(nn.getCntNonNei(v) + edC - sg.degreeDeep(deep, v) == k) {
                        // ok = false;
                        canAdd = false;
                        break;
                    }
                }
                if(!canAdd) break;
            }

            if(canAdd)
            for(ui j = 0, cntNN = nn.getCntNonNei(u); j < cntNN; j++) {
                ui v = nn.buffer[u*k+j];
                if(nn.getCntNonNei(v) + edC - sg.degreeDeep(deep, v) == k) {
                    canAdd = false;
                    break;
                }
            }

            if(canAdd) {
                ok = false;
                break;
            }
        }

        if(ok) {
#ifdef DDDEBUG
printf(" P+Cok\n");
#endif
            answer++;
        }

        return;
    }
#ifdef DDDEBUG
else printf(" no P+C ; minDPC %u\n", minDPC);
#endif

    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;
    
    ui countOfNonNeihborOfPv = nn.getCntNonNei(pv);
    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        if(g.connectHash(pv, v)) {
            // if(hasCommenNonNeighbor(pv, v)) cand[candSize++] = v;
            for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                if(!g.connectHash(v, nn.buffer[pv*k + i])) {
                    cand[candSize++] = v;
                    break;
                }
            }
        }
        else cand[candSize++] = v;
    }
#ifdef DDDEBUG
printf("cand:");
for(ui j = 0; j < candSize; j++) printf("%u ", cand[j]);
printf("\n");
#endif

    auto delDegreeInC = [&](ui v) {
#ifdef DDDEBUG
printf("del %u:", v);
for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) printf("%u ", C[i]);
printf("\n");
for(ui i = stX; i < edX; i++) if(g.connectHash(v, X[i])) printf("%u ", X[i]);
printf("\n");
for(ui i = 0; i < edP; i++) if(g.connectHash(v, P[i])) printf("%u ", P[i]);
printf("\n");
#endif
        // if(edC + edP + edX-stX < sg.degree2(v)) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) sg.deg[deep][C[i]]--;
            for(ui i = stX; i < edX; i++) if(g.connectHash(v, X[i])) sg.deg[deep][X[i]]--;
            for(ui i = 0; i < edP; i++) if(g.connectHash(v, P[i])) sg.deg[deep][P[i]]--;
        // }
//         else
//         for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
//             ui u = sg.pEdge[i];
// #ifdef DDDEBUG
// printf("s%u ", u);
// #endif
//             if(C.isIn(u, 0, edC) || X.isIn(u, stX, edX) || P.isIn(u, 0, edP)) {
// #ifdef DDDEBUG
// printf("%u ", u);
// #endif
//                 sg.deg[deep][u]--;
//             }
//         }
// #ifdef DDDEBUG
// printf("\n");
// #endif
    };
    auto addDegreeInC = [&](ui v) {
        // if(edC + edP + edX-stX < sg.degree2(v)) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) sg.deg[deep][C[i]]++;
            for(ui i = stX; i < edX; i++) if(g.connectHash(v, X[i])) sg.deg[deep][X[i]]++;
            for(ui i = 0; i < edP; i++) if(g.connectHash(v, P[i])) sg.deg[deep][P[i]]++;
        // }
        // else
        // for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
        //     ui v = sg.pEdge[i];
        //     if(C.isIn(v, 0, edC) || X.isIn(v, stX, edX) || P.isIn(v, 0, edP)) {
        //         sg.deg[deep][v]++;
        //     }
        // }
    };

    auto updateC = [&](ui v, LinearSet & T, ui ed) {
        ui newEnd = ed;
        ui cntNN = nn.getCntNonNei(v);

        for(ui j = 0; j < cntNN; j++) {
            ui w = nn.buffer[v*k + j];
            if(nn.getCntNonNei(w) == k) {
                ui tmp = 0;
                if(sg.degree2(w) <= newEnd) {
                    for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) {
                        if(T.isIn(sg.pEdge[l], 0, newEnd)) {
                            T.changeTo(sg.pEdge[l], tmp++);
                        }
                    }
                }
                else {
                    for(ui l = 0; l < newEnd; l++) {
                        if(g.connectHash(w, T[l])) {
                            T.changeTo(T[l], tmp++);
                        }
                    }
                }

                for(ui i = tmp; i < newEnd; i++) delDegreeInC(T[i]);

                newEnd = tmp;
            }
        }

        return newEnd;
    };

    auto updateX = [&](ui v, LinearSet & T, ui st, ui ed) {
        ui newSt = st;
        ui cntNN = nn.getCntNonNei(v);

        for(ui j = 0; j < cntNN; j++) {
            ui w = nn.buffer[v*k + j];;
// #ifdef DDDEBUG
// printf("upX:%u-%u:", v, w);
// #endif
            if(nn.getCntNonNei(w) == k) {
//                 if(sg.degree2(w) <= (ed - newSt)) {
//                     ui tmp = ed;
//                     for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) {
// printf("tt%u-%u %u\n", v, w, sg.pEdge[l]);
//                         if(T.isIn(sg.pEdge[l], newSt, ed)) {

// printf("tt%u-%u %u\n", v, w, sg.pEdge[l]);

//                             T.changeTo(sg.pEdge[l], --tmp);
//                         }
//                     }
// #ifdef DDDEBUG
// for(ui ii = newSt; ii < tmp; ii++)
// printf("%u-%d ", T[ii], (int)g.connectHash(w, T[ii]));
// printf("\n");
// #endif
//                     newSt = tmp;
//                 }
//                 else {
                    ui tmp = ed;
                    for(ui l = newSt; l < tmp; ) {
                        if(g.connectHash(w, T[l])) T.changeToByPos(l, --tmp);
                        else l++;
                    }
// #ifdef DDDEBUG
// for(ui ii = newSt; ii < tmp; ii++)
// printf("%u+%d ", T[ii], (int)g.connectHash(w, T[ii]));
// printf("\n");
// #endif
                    newSt = tmp;
                // }
            }
        }

        return newSt;
    };

    // for(ui i = 0; i < edC; i++) sg.pDeepIdx[deep+1][C[i]] = sg.pDeepIdx[deep][C[i]];
    // for(ui i = stX; i < edX; i++) sg.pDeepIdx[deep+1][X[i]] = sg.pDeepIdx[deep][X[i]];
    // for(ui i = 0; i < edP; i++) sg.pDeepIdx[deep+1][P[i]] = sg.pDeepIdx[deep][P[i]];
    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];

        //build next C 
        //build next X
        /**
         * if v's non-neighbor set W in P has count k, limit C to be W's neighbors
         * else:
         * delete the neighbors of v that has non-neibors > k in P
         * equal to
         * delete the neighbors of v that has neighbors <= |P|-k in P
        **/
#ifdef DDDEBUG
printf("    deep %u-%u, cand %u\n", deep, i, v);
#endif
#ifdef DDDEBUG
printf("deep %u\n", deep);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = stX; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
#endif
        if(nn.getCntNonNei(v) >= k) {
#ifdef DDDEBUG
printf("continue %u\n", v);
#endif
            C.changeTo(v, --edC);
            delDegreeInC(v);
            X.changeTo(v, edX++);
            continue;
        }

        //the prune technique of Dai
        for(ui j = 0; j < k; j++) bucket[j].clear();
        if(sg.degree2(v) < edC)
        for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
            ui u = sg.pEdge[j];
            if(C.isIn(u, 0, edC)) {
                ui cNei = nn.getCntNonNei(u);
                bucket[cNei].push_back(u);
            }
        }
        else {
            for(ui j = 0; j < edC; j++) {
                if(g.connectHash(v, C[j])) {
                    ui cNei = nn.getCntNonNei(C[j]);
                    bucket[cNei].push_back(C[j]);  
                }
            }
        }
        if(support.size() < edP) support.resize(edP * 2+1);
        ui s = 0;
        for(ui j = 0; j < edP; j++) {
            support[j] = k - nn.getCntNonNei(P[j]);
            s += k - nn.getCntNonNei(P[j]);
        }
        ui up = edP + k - nn.getCntNonNei(v)+ bucket[0].size();
        for(ui j = 1; j < k; j++) {
            for(auto u: bucket[j]) {
                ui minJ = 0, minS = k + 1;
                for(ui l = 0; l < edP; l++) {
                    if(!g.connectHash(u, P[l])) {
                        if(support[l] < minS) {
                            minJ = l;
                            minS = support[l];
                        }
                    }
                }
// assert(minS < k + 1);
                if(minS > 0 && minS < k) {
                    up++;
                    s -= j;
                    support[minJ]--;
                }

                if(s < j) break;
            }
            if(s < j) break;
        }
        if(up < q) {
#ifdef DDDEBUG
printf("up prune\n");
#endif
            C.changeTo(v, --edC);
            delDegreeInC(v);
            X.changeTo(v, edX++);
            continue;
        }

        P.changeTo(v, edP++);
        C.changeTo(v, --edC);
        delDegreeInC(v);

        //build tmpC, tmpStX and tmpEdX for the next level
        ui tmpC = edC, tmpStX = stX, tmpEdX = edX;
// ui tt = 0;
//         for(ui j = 0; j < edP; j++) if(!g.connectHash(v, P[j])) {
//             tt++,nn.addNonNei(P[j], v);
// int f = 0;
// for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
//     if(P[j] == nn.buffer[v*k + i]) f++;
// }
// assert(f == 1);
//         }
// assert(nn.getCntNonNei(v) == tt);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.addNonNei(nn.buffer[v*k + i], v);
        }
        nn.addNonNei(v, v);
        for(ui j = 0; j < edC; j++) if(!g.connectHash(v, C[j])) nn.addNonNei(C[j], v);
        for(ui j = stX; j < edX; j++) if(!g.connectHash(v, X[j])) nn.addNonNei(X[j], v);

        //if v in P has k non-neighors
        tmpC = updateC(v, C, edC);
        tmpStX = updateX(v, X, stX, edX);
        //setp 2
        //remove non-neighbor of v that have more than k non-neighbros
        for(ui j = 0; j < tmpC; ) {
            if(nn.getCntNonNei(C[j]) >= k) {
                delDegreeInC(C[j]);
                C.changeToByPos(j, --tmpC);
            }
            else j++;
        }
        for(ui j = tmpStX; j < edX; j++) {
            if(nn.getCntNonNei(X[j]) >= k) {
                X.changeToByPos(j, tmpStX++);
            }
        }
// if(v == pv) {
//     ui a = 0, b = 0, c = 0;
//     for(ui j = 0; j < tmpC; j++) {
//         if(g.connectHash(pv, C[j])) {
//             bool is = false;
//             for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
//                 if(g.connectHash(C[j], nn.buffer[pv*k + i])) {
//                     is = true;
//                     break;
//                 }
//             }
//             if(is) b++;
//             else a++;
//         }
//         else c++;
//     }
//     printf("%u, %u %u %u\n", deep, a, b, c);
// }
        bkPivot(deep + 1, tmpStX, edX, tmpC, edP);

        for(ui j = tmpC; j < edC; j++) addDegreeInC(C[j]);
        // for(ui j = 0; j < edP; j++) if(!g.connectHash(v, P[j])) nn.pop(P[j]);
        nn.pop(v);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.pop(nn.buffer[v*k + i]);
        }
        for(ui j = 0; j < edC; j++) if(!g.connectHash(v, C[j])) nn.pop(C[j]);
        for(ui j = stX; j < edX; j++) if(!g.connectHash(v, X[j])) nn.pop(X[j]);
        
        P.changeTo(v, --edP);
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
