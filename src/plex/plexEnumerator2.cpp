/**
 * choose the pivot in C and X 
*/

#include "plexEnumerator.h"
// #define DDDEBUG


ui plexEnumerator::run() {
    printf("plexEnumerator2\n");
    std::vector<ui> toDelete(g.maxD);
    toDelete.clear();
#ifdef DDDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<' '<<cnt<<std::endl;
#endif
        //P is empty
        //X is 2-hop neighbors < u
        //C is 2-hop neighbors > u
        if(g.pIdx[u + 1] - g.pIdx2[u] + k < q) continue;
        edC = 0;
        stX = edX = 0;
        neighborsInP.clear();
        neighborsInP.insert(u, 0);
        edP = 0;
        P.changeTo(u, edP++);

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
                        if(g.connect(v, w)) {
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
                    if(g.connect(v, w)) {
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
        for(ui j = 0; j < edC; j++) neighborsInP.insert(C[j], 1);
        for(ui j = 0; j < edX; j++) neighborsInP.insert(X[j], 1);
#ifdef DDDEBUG
printf("neighbors in P:\n");
for(ui i = 0; i < edC; i++) printf("%u ", neighborsInP[C[i]]); printf("\n");
for(ui i = 0; i < edX; i++) printf("%u ", neighborsInP[X[i]]); printf("\n");
#endif
        if(k > 1)
        for(ui j = 0, ed = edC; j < ed; j++) {
            ui v = C[j];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connect(u, w)) continue;
                if(C.isIn(w, 0, edC) || X.isIn(w, 0, edX) || w == u) continue; 
                ui d = 0;
                
                for(ui j = 0; j < edC; j++) {
                    ui x = C[j]; 
                    if(g.connect(w, x)) {
                        d++;

                        if(d + 2 * k >= q + 2) {
                            if(w > u) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in C\n", w, v);
#endif
                                C.changeTo(w, edC++);
                                // assert(g.connect(u, w) == false);
                                neighborsInP.insert(w, 0);
                            }
                            else if(w < u) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in X\n", w, v);
#endif
                                X.changeTo(w, edX++);
                                // assert(g.connect(u, w) == false);
                                neighborsInP.insert(w, 0);
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
#ifdef DDDEBUG
printf("%u:", v);
#endif
            for(ui i = g.pIdx[v]; i < g.pIdx[v + 1]; i++) {

                ui w = g.pEdge[i];

                if(C.isIn(w, 0, edC) || X.isIn(w, 0, edX) || w == u) {
                    sg.pEdge[sg.pIdx2[v]++] = w;
#ifdef DDDEBUG
printf("%u ", w);
#endif
                }
            }
#ifdef DDDEBUG
printf("\n");
#endif
        };
#ifdef DDDEBUG
printf("sg:\n");
#endif
        for(ui i = 0; i < edC; i++) buildV(C[i]);
        for(ui i = 0; i < edX; i++) buildV(X[i]);
        
        buildV(u);

#ifdef DDDEBUG
printf("neighbors in P:\n");
for(ui i = 0; i < edC; i++) printf("%u ", neighborsInP[C[i]]); printf("\n");
for(ui i = 0; i < edX; i++) printf("%u ", neighborsInP[X[i]]); printf("\n");
#endif

        bkPivot(0);
    }

    return cnt;
}

void plexEnumerator::bkPivot(ui deep) {
#ifdef DDDEBUG
printf("deep %u\n", deep);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("X:"); for(ui i = stX; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
#endif
    if(edC == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        if(stX == edX && edP >= q) {
            for(ui i = 0; i < edP; i++) {
                if(neighborsInP[P[i]] == 0) return; 
            }
#ifdef DDDEBUG
printf("    P is ans!\n");
#endif
            cnt++;
        }
        return;
    }

    if(edP + edC < q) return;

    //is P+C a k-plex
    ui minDPC = g.n;
    //find pivot in C and X with max degree in C
    ui maxD = 0, pv = C[0];
    bool inC = true;

    for(ui i = 0; i < edC; i++) {
        ui v = C[i];
        ui d = 0;

        if(sg.degree2(v) <= edC * logBasedOn2[sg.degree2(v)]) {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                if(C.isIn(sg.pEdge[i], 0, edC)) {
                    d++;
                }
            }
// printf("find d in C1:%u %u\n", v, d);
        }
        else {
            for(ui j = 0; j < edC; j++) {
                ui w = C[j];
                if(sg.connect2(v, w)) d++;
            }
// printf("find d in C2:%u %u\n", v, d);
        }
        
        minDPC = std::min(minDPC, d + neighborsInP[v]);
        if(d > maxD) {
            maxD = d;
            pv = v;
        }
    }
    for(ui i = stX; i < edX; i++) {
        ui v = X[i];
        ui d = 0;

        if(sg.degree2(v) <= edC * logBasedOn2[sg.degree2(v)]) {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                if(C.isIn(sg.pEdge[i], 0, edC)) {
                    d++;
                }
            }
// printf("find d in C1:%u %u\n", v, d);
        }
        else {
            for(ui j = 0; j < edC; j++) {
                ui w = C[j];
                if(sg.connect2(v, w)) d++;
            }
// printf("find d in C2:%u %u\n", v, d);
        }
        
        if(d > maxD) {
            maxD = d;
            pv = v;
            inC = false;
        }
    }

    //compute minDPC
    for(ui i = 0; i < edP; i++) {
        ui v = P[i];
        ui d = 0;

        if(sg.degree2(v) <= edC * logBasedOn2[sg.degree2(v)]) {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                if(C.isIn(sg.pEdge[i], 0, edC)) {
                    d++;
                }
            }
// printf("find d in C1:%u %u\n", v, d);
        }
        else {
            for(ui j = 0; j < edC; j++) {
                ui w = C[j];
                if(sg.connect2(v, w)) d++;
            }
// printf("find d in C2:%u %u\n", v, d);
        }
        
        minDPC = std::min(minDPC, d + neighborsInP[v]);
    }

#ifdef DDDEBUG
printf("minDPC: %u\n", minDPC);
printf("maxD %u, pv %u\n", maxD, pv);
#endif
    if(minDPC + k >= edP + edC) {//P+C is a k-plex, check X
        bool ok = true;

        for(ui i = 0; i < edC; i++) {
            
            // C.changeTo(v, --edC);
            // for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) neighborsInP[sg.pEdge[j]]++;
            ui v = C[i];
            P.changeTo(v, edP++);
            for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
                neighborsInP.increment(sg.pEdge[j]);
            }
        }

        //delete x 
        for(ui i = 0; i < edP; i++) {
            ui v = P[i];
            
            if(edP - neighborsInP[v] == k) {
                ui tmp = edX;
                for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
                    if(X.isIn(sg.pEdge[j], stX, edX)) {
                        X.changeTo(sg.pEdge[j], --tmp);
                    }
                }
                stX = tmp;
            }
        }
#ifdef DDDEBUG
printf("X:"); for(ui i = stX; i < edX; i++) printf("%u ", X[i]); printf("\n");
fflush(stdout);
#endif
        for(ui i = stX; i < edX; i++) {
            ui x = X[i];
            //if x+P+C is a k-plex
            ui nonNeighbor = edP - neighborsInP[x];
            if(nonNeighbor >= k) continue;
#ifdef DDDEBUG
printf("%u in X has %u nonNeighbors\n", x, nonNeighbor);
#endif
            if(nonNeighbor <= k - 1) {
#ifdef DDDEBUG
printf("    P+C is not ans!\n");
#endif
                ok = false;
                break;
            }
        }

        if(ok) {
#ifdef DDDEBUG
printf("    P+C is ans\n");
#endif
            cnt++;
        }

        for(ui i = 0; i < edC; i++) {
            // P.changeTo(v, edP++);
            // C.changeTo(v, --edC);
            // for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) neighborsInP[sg.pEdge[j]]++;
            ui v = C[i];
            P.changeTo(v, --edP);
            for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
                neighborsInP.decrement(sg.pEdge[j]);
            }
        }

        return;
    }

    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;
    
    ui countOfNonNeihborOfPv = 0;
    for(ui i = 0; i < edP; i++) {
        if(!sg.connect2(pv, P[i])) {
            nonNeighborsOfPv[countOfNonNeihborOfPv++] = P[i];
        }
    }
    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        if(sg.connect2(pv, v)) {
            // if(hasCommenNonNeighbor(pv, v)) cand[candSize++] = v;
            for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                if(!sg.connect2(v, nonNeighborsOfPv[i])) {
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
    // for(ui i = 0; i < edP; i++) {
    //     if(edP - neighborsInP[P[i]] == k - 1) {
            
    //     }
    // }

    auto updateC = [&](ui v, ui neiOfv, LinearSet & T, ui ed) {
        ui newEnd = ed;
#ifdef DDDEBUG
printf("updateC %u:\n", v);
#endif
        for(ui j = neiOfv; j < edP; j++) {
            ui w = P[j];
#ifdef DDDEBUG
printf("w:%u, edP:%u, neighborsInP[w]:%u, newEnd %u\n", 
    w, edP, neighborsInP[w], newEnd);
#endif
            if(edP - neighborsInP[w] == k) {
                ui tmp = 0;
                if(sg.degree2(w) <= newEnd * logBasedOn2[sg.degree2(w)]) {
#ifdef DDDEBUG
printf("w1:%u, deg %u:", w, sg.pIdx2[w] - sg.pIdx[w]);
for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) printf("%u ", sg.pEdge[l]);printf("\n");
#endif
                    for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) {
#ifdef DDDEBUG
printf("w:%u, %u\n", w, sg.pEdge[l]);
#endif
                        if(T.isIn(sg.pEdge[l], 0, newEnd)) {
                            T.changeTo(sg.pEdge[l], tmp++);
#ifdef DDDEBUG
printf("inter %u \n", sg.pEdge[l]);
#endif
                        }
#ifdef DDDEBUG
else
printf("not1 inter %u \n", sg.pEdge[l]);
#endif
                    }
                }
                else {
#ifdef DDDEBUG
printf("w2:%u, edP:%u, neighborsInP[w]:%u, newEnd %u\n", 
    w, edP, neighborsInP[w], newEnd);
#endif
                    for(ui l = 0; l < newEnd; l++) {
#ifdef DDDEBUG
printf("w:%u, %u\n", w, T[l]);
#endif
                        if(sg.connect2(w, T[l])) {
                            T.changeTo(T[l], tmp++);
#ifdef DDDEBUG
printf("inter %u \n", T[l]);
#endif
                        }
#ifdef DDDEBUG
else
printf("not2 inter %u \n", T[l]);
#endif
                    }
                }

                newEnd = tmp;
            }

        }

        return newEnd;
    };

    auto updateX = [&](ui v, ui neiOfv, LinearSet & T, ui st, ui ed) {
        ui newSt = st;
#ifdef DDDEBUG
printf("updateX %u:\n", v);
#endif
        for(ui j = neiOfv; j < edP; j++) {
            ui w = P[j];
#ifdef DDDEBUG
printf("w:%u, edP:%u, neighborsInP[w]:%u\n", w, edP, neighborsInP[w]);
#endif
            if(edP - neighborsInP[w] == k) {
                if(sg.degree2(w) <= (ed - newSt) * logBasedOn2[sg.degree2(w)]) {
                    ui tmp = ed;
                    for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) {
                        if(T.isIn(sg.pEdge[l], newSt, ed)) T.changeTo(sg.pEdge[l], --tmp);
                    }
                    newSt = tmp;
                }
                else {
                    ui tmp = ed;
                    for(ui l = newSt; l < tmp; ) {
                        if(sg.connect2(w, T[l])) T.changeTo(T[l], --tmp);
                        else l++;
                    }
                    newSt = tmp;
                }
            }
        }

        return newSt;
    };

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
        if(edP - neighborsInP[v] >= k) {
#ifdef DDDEBUG
printf("continue %u\n", v);
#endif
            continue;
        }

        P.changeTo(v, edP++);
        C.changeTo(v, --edC);
        // for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) neighborsInP[sg.pEdge[j]]++;
        for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
            neighborsInP.increment(sg.pEdge[j]);
        }
        
        ui neiOfv = 0;
        if(sg.degree2(v) <= edP * logBasedOn2[sg.degree2(v)]) {
            for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
                if(P.isIn(sg.pEdge[j], 0, edP)) P.changeTo(sg.pEdge[j], neiOfv++);
            }
        }
        else {
            for(ui j = 0; j < edP; j++) {
                if(sg.connect2(v, P[j])) {
                    P.changeTo(P[j], neiOfv++);
                }
            }
        }
#ifdef DDDEBUG

printf("deep%u: %u nei/non-nei in P:", deep, v);
for(ui j = 0; j < neiOfv; j++) printf("%u ", P[j]); printf(" / ");
for(ui j = neiOfv; j < edP; j++) printf("%u ", P[j]); printf("\n");
#endif
        ui tmpC = edC, tmpStX = stX, tmpEdX = edX;
        edC = updateC(v, neiOfv, C, edC);
        stX = updateX(v, neiOfv, X, stX, edX);

        //setp 2
        //remove non-neighbor of v that have more than k non-neighbros
        for(ui j = 0; j < edC; ) {
#ifdef DDDEBUG
printf("step2: j %u, C[j] %u:", j, C[j]);
#endif
            // if(!sg.connect2(v, C[j])) {
#ifdef DDDEBUG

printf("connected ");
#endif
                if(edP - neighborsInP[C[j]] >= k) {
#ifdef DDDEBUG

printf("edP %u, %u, %u", edP, C[j], neighborsInP[C[j]]);
#endif
                    C.changeToByPos(j, --edC);
                }
                else j++;
            // }
            // else j++;
#ifdef DDDEBUG

printf("\n");
#endif
        }

        for(ui j = stX; j < edX; j++) {
#ifdef DDDEBUG

printf("step2: j%u, X[j]%u:", X[j]);
#endif
           // if(!sg.connect2(v, X[j])) {
// #ifdef DDDEBUG

// printf("connected ");
// #endif
                if(edP - neighborsInP[X[j]] >= k) {
#ifdef DDDEBUG

printf("edP %u, %u, %u", edP, X[j], neighborsInP[X[j]]);
#endif
                    X.changeToByPos(j, stX++);
                }
          //  }
          //  else j++;
#ifdef DDDEBUG

printf("\n");
#endif
        }

        bkPivot(deep + 1);
        edC = tmpC, stX = tmpStX, edX = tmpEdX; 
        
        P.changeTo(v, --edP);
        X.changeTo(v, edX++);
        // for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) neighborsInP[sg.pEdge[j]]--;
        for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
            neighborsInP.decrement(sg.pEdge[j]);
        }
    }

    if(candSize > 0) {
        for(ui i = candSize - 1; i > 0; i--) {
            ui v = cand[i];
            X.changeTo(v, --edX);
        }
        X.changeTo(cand[0], --edX);
    }
    
}

