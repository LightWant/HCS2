#include "plexCounting.h"
#include <iostream>

// #define BASELINE
#ifdef BASELINE
#define BASELINEK 3
#define BASELINEQ 10
#endif

// #define DDDEBUG

std::vector<ull> plexCounting::run() {
#ifdef DDDEBUG
g.print();
#endif
    printf("plexCounting.cpp::run\n");
    
    std::vector<ui> toDelete(g.maxD);
    toDelete.clear();

    g.initHash();
    printf("init Hash\n");

    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<' '<<answers[q]<<std::endl;
#endif
        //P is empty
        //X is 2-hop neighbors < u
        //C is 2-hop neighbors > u
        if(g.pIdx[u + 1] - g.pIdx2[u] + k < q) continue;
        ui edC = 0;
        ui edP = 0;

        P.changeTo(u, edP++);
        nn.addNonNei(u, u);

        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) C.changeTo(g.pEdge[i], edC++);
#ifdef DDDEBUG
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("edC %u\n", edC);
#endif
        //delete C
        while(true) {
            for(ui i = 0; i < edC; i++) {
                ui v = C[i];
                ui d = 0;
                if(g.degree(v) <= edC) {
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
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
#endif
        //two hop neighbors
        if(k > 1)
        for(ui j = 0, ed = edC; j < ed; j++) {
            ui v = C[j];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(w < u || g.connectHash(u, w)) continue;
                if(C.isIn(w, 0, edC) || w == u) continue; 
                ui d = 0;
                
                for(ui j = 0; j < ed; j++) {
                    ui x = C[j]; 
                    if(g.connectHash(w, x)) {
                        d++;

                        if(d + 2 * k >= q + 2) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in C\n", w, v);
#endif
                            C.changeTo(w, edC++);
                            nn.addNonNei(w, u);
                            break;
                        }
                    }
                }
            }
        }

        
#ifdef DDDEBUG
printf("After Add 2-hop neighbors\n");
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
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
        for(ui i = 0; i < edC; i++) sg.deg[0][C[i]] = sg.pIdx2[C[i]] - sg.pIdx[C[i]];
        buildV(u);
        sg.deg[0][u] = sg.pIdx2[u] - sg.pIdx[u];
        for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
            ui w = g.pEdge[i];
            if(C.isIn(w, 0, edC)) {
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
#endif

        bkPivot(1, edC, edP, 0, 1);

        for(ui i = 0; i < edC; i++) nn.clear(C[i]);
        for(ui i = 0; i < edP; i++) nn.clear(P[i]);
    }
#ifdef BASELINE
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != BASELINEQ) return false;

    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
        ui d = 0;
        for(ui v = 0; v < g.n; v++) if((1<<v) & x){
            if(g.connectHash(u, v)) d++;
        }
        if(d + BASELINEK < BASELINEQ) return false; 
    }
    return true;
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

void plexCounting::bkPivot(ui deep, ui edC, ui edP, ui p, ui h) {
#ifdef DDDEBUG
printf("deep %u, p %u, h %u\n", deep, p, h);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
#endif
    if(edC == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        if(p + h >= q) {
#ifdef DDDEBUG
printf("    P is ans!\n");
#endif
            for(ui i = std::max(q, h); i <= p + h && i < maxSize; i++) {
                answers[i] += CN[p][i - h];
            }
        }
        return;
    }

    if(p + h + edC < q) return;

    //init deg from the pre level
    if(sg.deg[deep].size() == 0) sg.deg[deep].resize(g.n);
    for(ui i = 0; i < edC; i++) sg.deg[deep][C[i]]=sg.deg[deep-1][C[i]];
    for(ui i = 0; i < edP; i++) sg.deg[deep][P[i]]=sg.deg[deep-1][P[i]];
#ifdef DDDEBUG
printf("CdeginC:");
for(ui i = 0; i < edC; i++) printf("%u-%u ", C[i], sg.deg[deep][C[i]]);
printf("\nPdeginC:");
for(ui i = 0; i < edP; i++) printf("%u-%u ", P[i], sg.deg[deep][P[i]]);
printf("\n");
#endif

    //is P+C a k-plex, minDPC is the min degree in P+C
    ui minDPC = sg.deg[deep][C[0]] + edP - nn.getCntNonNei(C[0]);
    //find pivot in C with max degree in C
    ui maxD = sg.deg[deep][C[0]], pv = C[0];
    for(ui i = 1; i < edC; i++) {
        ui v = C[i];
        ui d = sg.degreeDeep(deep, v);
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
        if(d > maxD) {
            maxD = d;
            pv = v;
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

    if(minDPC + k >= edP + edC) {//P+C is a k-plex,
#ifdef DDDEBUG
printf(" P+C!; minDPC %u\n", minDPC);
#endif
        bkPivot(deep + 1, 0, 0, p + edC, h);
        return;
    }

    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;
    
    bool isP = true;
    ui countOfNonNeihborOfPv = nn.getCntNonNei(pv);
    ui C1Size = 0, C2Size = 0;
    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        if(v == pv) continue;
        if(g.connectHash(pv, v)) {
            // if(hasCommenNonNeighbor(pv, v)) cand[candSize++] = v;
            bool isInC1 = true;
            for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                ui w = nn.buffer[pv*k + i];
                if(!g.connectHash(v, w)) {
                    cand[candSize++] = v;
                    C2Size++;
                    isInC1 = false;
                    break;
                }
            }
            if(isInC1) C.changeToByPos(j, C1Size++);
            // if(isInC1) C1Size++;
        }
        else {
            // if(canAdd(v)) C2Size++;
            if(countOfNonNeihborOfPv+1 < k && nn.getCntNonNei(v)+1 < k) {
                bool canAdd = true;
                for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                    ui w = nn.buffer[pv*k + i];
                    if(nn.getCntNonNei(w)+1 == k && !g.connectHash(v, w)) {
                        canAdd = false;
                    }
                }
                if(canAdd) C2Size++;
            }
            
            cand[candSize++] = v;
        }
    }
#ifdef DDDEBUG
printf("pv: %u", pv);
printf("cand:");
for(ui i = 0; i < candSize; i++) printf("%u ", cand[i]);printf("\n");
#endif

    auto delDegreeInC = [&](ui v) {
#ifdef DDDEBUG
printf("del: %u:", v);
#endif
        if(edC + edP < sg.degree2(v)) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) sg.deg[deep][C[i]]--;
            for(ui i = 0; i < edP; i++) if(g.connectHash(v, P[i])) sg.deg[deep][P[i]]--;
        }
        else
        for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
            ui u = sg.pEdge[i];

            if(C.isIn(u, 0, edC) || P.isIn(u, 0, edP)) {
#ifdef DDDEBUG
printf("%u ", u);
#endif
                sg.deg[deep][u]--;
            }
        }
#ifdef DDDEBUG
printf("\n");
#endif
    };
    auto addDegreeInC = [&](ui v) {
        if(edC + edP < sg.degree2(v)) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) sg.deg[deep][C[i]]++;
            for(ui i = 0; i < edP; i++) if(g.connectHash(v, P[i])) sg.deg[deep][P[i]]++;
        }
        else
        for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
            ui v = sg.pEdge[i];
            if(C.isIn(v, 0, edC) || P.isIn(v, 0, edP)) {
                sg.deg[deep][v]++;
            }
        }
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

    auto solvePre = [&](ui v) {
        P.changeTo(v, edP++);
        C.changeTo(v, --edC);
        delDegreeInC(v);

        // for(ui j = 0; j < edP; j++) if(!g.connectHash(v, P[j])) nn.addNonNei(P[j], v);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.addNonNei(nn.buffer[v*k + i], v);
        }
        nn.addNonNei(v, v);
        for(ui j = 0; j < edC; j++) if(!g.connectHash(v, C[j])) nn.addNonNei(C[j], v);

        //build tmpC for the next level
        //if v in P has k non-neighors
        ui tmpC = updateC(v, C, edC);
        //setp 2
        //remove non-neighbor of v that have more than k non-neighbros
        for(ui j = 0; j < tmpC; ) {
            if(nn.getCntNonNei(C[j]) >= k) {
                delDegreeInC(C[j]);
                C.changeToByPos(j, --tmpC);
            }
            else j++;
        }

        return tmpC;
    };

    auto solveBack = [&](ui v, ui newEdC) {
        for(ui j = newEdC; j < edC; j++) addDegreeInC(C[j]);
        // for(ui j = 0; j < edP; j++) if(!g.connectHash(pv, P[j])) nn.pop(P[j]);
        nn.pop(v);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.pop(nn.buffer[v*k + i]);
        }
        for(ui j = 0; j < edC; j++) if(!g.connectHash(v, C[j])) nn.pop(C[j]);
        P.changeTo(pv, --edP);
    };

    if(C2Size == 0) {
#ifdef DDDEBUG
printf("    deep %u-%u, C2=0\n", deep, pv);
#endif   
        ui newEdC = solvePre(pv);
        bkPivot(deep + 1, newEdC, edP, p + 1, h);
        solveBack(pv, newEdC);
    }
    else  {
#ifdef DDDEBUG
printf("    deep %u-%u, C2>0\n", deep, pv);
#endif  
        P.changeTo(pv, edP++);
        for(ui i = C1Size; i < edC; i++) delDegreeInC(C[i]);
        bkPivot(deep + 1, C1Size, edP, p, h);
        for(ui i = C1Size; i < edC; i++) addDegreeInC(C[i]);
        P.changeTo(pv, --edP);

        ui newEdC = solvePre(pv);
        bkPivot(deep + 1, newEdC, edP, p, h + 1);
        solveBack(pv, newEdC);
    }

    auto upPrune = [&](ui v) {
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
// printf("cNei %u %u\n", cNei, C[j]);fflush(stdout);
// assert(cNei < k); 
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
        ui up = edP + std::min(edC - sg.deg[deep][v], k - nn.getCntNonNei(v)) + bucket[0].size();
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
        
        return up < q;
    };

    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];
#ifdef DDDEBUG
printf("    deep %u-%u, cand %u\n", deep, i, v);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
fflush(stdout);
#endif   
        if(nn.getCntNonNei(v) >= k) {
#ifdef DDDEBUG
printf("continue %u\n", v);
#endif
            C.changeTo(v, --edC);
            delDegreeInC(v);
            continue;
        }

        if(upPrune(v)) {
            C.changeTo(v, --edC);
            delDegreeInC(v);
            continue;
        }

        ui newEdC = solvePre(v);
        bkPivot(deep + 1, newEdC, edP, p, h + 1);
        solveBack(v, newEdC);
    }
}
