#include "d2k.h"
// #define DDDEBUG

ui d2k::run() {
#ifdef DDDEBUG
g.print();
#endif
    std::vector<ui> toDelete(g.maxD);
    toDelete.clear();

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
        //two hop neighbors
        if(k > 1)
        for(ui j = 0, ed = edC; j < ed; j++) {
            ui v = C[j];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(g.connect(u, w)) continue;
                if(C.isIn(w, 0, edC) || X.isIn(w, 0, edX) || w == u) continue; 
                ui d = 0;
                
                for(ui j = 0; j < ed; j++) {
                    ui x = C[j]; 
                    if(g.connect(w, x)) {
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

        bkPivot(0, stX, edX, edC, edP);

        for(ui i = 0; i < edC; i++) nn.clear(C[i]);
        for(ui i = 0; i < edP; i++) nn.clear(P[i]);
        for(ui i = stX; i < edX; i++) nn.clear(X[i]);
    }

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
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
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
        }
        else {
            for(ui j = 0; j < edC; j++) {
                ui w = C[j];
                if(sg.connect2(v, w)) d++;
            }
        }
        
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
        ui d = 0;

        if(sg.degree2(v) <= edC * logBasedOn2[sg.degree2(v)]) {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                if(C.isIn(sg.pEdge[i], 0, edC)) {
                    d++;
                }
            }
        }
        else {
            for(ui j = 0; j < edC; j++) {
                ui w = C[j];
                if(sg.connect2(v, w)) d++;
            }
        }
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
        if(minDPC + k < edP + edC) break;
    }

    if(minDPC + k >= edP + edC) {//P+C is a k-plex, check X
#ifdef DDDEBUG
printf(" P+C!; minDPC %u\n", minDPC);
#endif
        bool ok = true;

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
            
            support[i] = d;
// printf("supportC:%u %u\n", v, d);
        }
        for(ui i = 0; i < edP; i++) {
            ui v = P[i];
            ui d = 0;

            if(sg.degree2(v) <= edC * logBasedOn2[sg.degree2(v)]) {
                for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                    if(C.isIn(sg.pEdge[i], 0, edC)) {
                        d++;
                    }
                }
            }
            else {
                for(ui j = 0; j < edC; j++) {
                    ui w = C[j];
                    if(sg.connect2(v, w)) d++;
                }
            }
            
            support[edC + i] = d;
// printf("supportD:%u %u\n", v, d);
        }

        for(ui i = stX; i < edX; i++) {
            ui u = X[i];
            bool canAdd = true;

            //check if P+C+u is a plex
            ui nnU = nn.getCntNonNei(u) + 1;
// printf("X %u, nnU %u\n", u, nnU);
            if(nnU <= k)
            for(ui j = 0; j < edC; j++) {
                ui v = C[j];
// printf("C %u\n", v);
                if(!sg.connect2(u, v)) {

                    nnU++;
// printf("nonNei, nnU %u, %u %u %u\n", nnU, nn.getCntNonNei(v), edC, support[j]);
                    if(nnU > k) {
                        canAdd = false;
                        break;
                    }
                    if(nn.getCntNonNei(v) + edC - support[j] == k) {
                        // ok = false;
                        canAdd = false;
                        break;
                    }
                }
                if(!canAdd) break;
            }

            if(nnU <= k && canAdd)
            for(ui j = 0, cntNN = nn.getCntNonNei(u); j < cntNN; j++) {
                ui v = nn.buffer[u*k+j];
                ui idxV = P.idx(v);
                if(nn.getCntNonNei(v) + edC - support[edC+idxV] == k) {
                    canAdd = false;
                    break;
                }
            }

            if(nnU > k) canAdd = false;

            if(canAdd) ok = false;

            if(!ok) break;
        }

        if(ok) {
#ifdef DDDEBUG
printf(" P+Cok\n");
#endif
            answer++;
        }

        return;
    }

    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;
    
    ui countOfNonNeihborOfPv = nn.getCntNonNei(pv);
    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        if(sg.connect2(pv, v)) {
            // if(hasCommenNonNeighbor(pv, v)) cand[candSize++] = v;
            for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                if(!sg.connect2(v, nn.buffer[pv*k + i])) {
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

    auto updateC = [&](ui v, LinearSet & T, ui ed) {
        ui newEnd = ed;
        ui cntNN = nn.getCntNonNei(v);

        for(ui j = 0; j < cntNN; j++) {
            ui w = nn.buffer[v*k + j];
            if(nn.getCntNonNei(w) == k) {
                ui tmp = 0;
                if(sg.degree2(w) <= newEnd * logBasedOn2[sg.degree2(w)]) {
                    for(ui l = sg.pIdx[w]; l < sg.pIdx2[w]; l++) {
                        if(T.isIn(sg.pEdge[l], 0, newEnd)) {
                            T.changeTo(sg.pEdge[l], tmp++);
                        }
                    }
                }
                else {
                    for(ui l = 0; l < newEnd; l++) {
                        if(sg.connect2(w, T[l])) {
                            T.changeTo(T[l], tmp++);
                        }
                    }
                }

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

            if(nn.getCntNonNei(w) == k) {
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
        if(nn.getCntNonNei(v) >= k) {
#ifdef DDDEBUG
printf("continue %u\n", v);
#endif
            C.changeTo(v, --edC);
            X.changeTo(v, edX++);
            continue;
        }

        //the prune technique of Dai
        for(ui j = 0; j < k; j++) bucket[j].clear();
        for(ui j = sg.pIdx[v]; j < sg.pIdx2[v]; j++) {
            ui u = sg.pEdge[j];
            if(C.isIn(u, 0, edC)) {
                ui cNei = nn.getCntNonNei(u);
                bucket[cNei].push_back(u);
            }
        }
        if(support.size() < edP) support.resize(edP * 2+1);
        ui s = 0;
        for(ui j = 0; j < edP; j++) {
            support[j] = k - nn.getCntNonNei(P[j]);
            s += k - nn.getCntNonNei(P[j]);
        }
        ui up = edP + k - nn.getCntNonNei(v) + bucket[0].size();
        for(ui j = 1; j < k; j++) {
            for(auto u: bucket[j]) {
                ui minJ = 0, minS = k + 1;
                for(ui l = 0; l < edP; l++) {
                    if(!sg.connect2(u, P[l])) {
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
            C.changeTo(v, --edC);
            X.changeTo(v, edX++);
            continue;
        }

        P.changeTo(v, edP++);
        C.changeTo(v, --edC);
        
        //build tmpC, tmpStX and tmpEdX for the next level
        ui tmpC = edC, tmpStX = stX, tmpEdX = edX;
        for(ui j = 0; j < edP; j++) if(!sg.connect2(v, P[j])) nn.addNonNei(P[j], v);
        for(ui j = 0; j < edC; j++) if(!sg.connect2(v, C[j])) nn.addNonNei(C[j], v);
        for(ui j = stX; j < edX; j++) if(!sg.connect2(v, X[j])) nn.addNonNei(X[j], v);

        
        tmpC = updateC(v, C, edC);
        tmpStX = updateX(v, X, stX, edX);

        //setp 2
        //remove non-neighbor of v that have more than k non-neighbros
        for(ui j = 0; j < tmpC; ) {
            if(nn.getCntNonNei(C[j]) >= k) {
                C.changeToByPos(j, --tmpC);
            }
            else j++;
        }
        for(ui j = tmpStX; j < edX; j++) {
            if(nn.getCntNonNei(X[j]) >= k) {
                X.changeToByPos(j, tmpStX++);
            }
        }

        bkPivot(deep + 1, tmpStX, edX, tmpC, edP);

        for(ui j = 0; j < edP; j++) if(!sg.connect2(v, P[j])) nn.pop(P[j]);
        for(ui j = 0; j < edC; j++) if(!sg.connect2(v, C[j])) nn.pop(C[j]);
        for(ui j = stX; j < edX; j++) if(!sg.connect2(v, X[j])) nn.pop(X[j]);
        
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
