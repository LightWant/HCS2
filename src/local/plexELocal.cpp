#include "plexLocal.h"

// #define BASELINE
// #define DDDEBUG
#ifdef DDDEBUG
#include <iostream>
#endif

void plexLocal::run() {
    printf("plexELocalv10::run\n");
#ifdef DDDEBUG
g.print();
#endif
    std::vector<ui> toDelete(std::min(g.n, g.maxD * k));
    toDelete.clear();

    g.initHash();

    nadj.resize(g.n);

    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<std::endl;
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
                            C.changeTo(w, edC++);
                            nn.addNonNei(w, u);
                            break;
                        }
                    }
                }
            }
        }

        //build sub-graph g
        auto buildV = [&](ui v) {
            nadj[v].clear();
            sg.pIdx[v] = sg.pIdx2[v] = g.pIdx[v];
            if(g.connectHash(v, u)) sg.pEdge[sg.pIdx2[v]++] = u;
            else nadj[v].push_back(u);
            for(ui i = 0; i < edC; i++) {
                if(g.connectHash(v, C[i])) sg.pEdge[sg.pIdx2[v]++] = C[i];
                else nadj[v].push_back(C[i]);
            }
        };
        for(ui i = 0; i < edC; i++) buildV(C[i]);
        
        buildV(u);

H[0] = u;
        bkPivot(1, edC, edP, 0, 1);

        for(ui i = 0; i < edC; i++) nn.clear(C[i]);
        for(ui i = 0; i < edP; i++) nn.clear(P[i]);
    }

    // for(ui i = 0; i < g.m/2; i++) {
    //     printf("%.0f\n", answers[reEegeId[i]]);
    // }

#ifdef BASELINE
std::vector<ui> ec(g.m/2);
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
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != q) return false;

    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
        ui d = 0;
        for(ui v = 0; v < g.n; v++) if((1<<v) & x){
            if(g.connectHash(u, v)) d++;
        }
        if(d + k < q) return false; 
    }
    return true;
};

ui cnt = 0;
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        cnt++;
        print(i);
// printf("xxx%u\n", ec[2]);
    }
}
printf("cnt:%u\n", cnt);
for(ui e = 0; e < g.m/2; e++) printf("e_cnt %u, %u %.0f\n", e, ec[e], answers[e]);
for(ui u = 0; u < g.n; u++) {
    for(ui e = pOIdx[u]; e < pOIdx[u+1]; e++) {
        printf("e%u, %u %u\n", e, u, pOEdge[e]);
    }
}
#endif
}

//Pn is the set of pivot, P is the enumerated nodes.
void plexLocal::bkPivot(ui deep, ui edC, ui edP, ui p, ui h) {
#ifdef DDDEBUG
printf("deep %u, p %u, h %u\n", deep, p, h);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("Pn:"); for(ui i = 0; i < p; i++) printf("%u ", Pn[i]); printf("\n");
printf("H:"); for(ui i = 0; i < h; i++) printf("%u ", H[i]); printf("\n");
printf("PP:"); for(auto e:PP) printf("%u ", e); printf("\n");
printf("\n");
fflush(stdout);
#endif
    auto updateAns = [&]() {
#ifdef DDDEBUG
printf("updateAns\n");
#endif
        if(q >= h+2)
        for(auto i : PP) {
            answers[i] += CN[p-2][q-h-2];
#ifdef DDDEBUG
printf("ans, %u, +%.0f pp\n", i, CN[p-2][q-h-2]);
#endif
        }
        if(q >= h+1)
        for(auto i : PH) {
#ifdef DDDEBUG
printf("ans, %u, +%.0f ph\n", i, CN[p-1][q-h-1]);
#endif
            answers[i] += CN[p-1][q-h-1];
        }
        for(auto i : HH) {
#ifdef DDDEBUG
printf("ans, %u, +%.0f hh\n", i, CN[p][q-h]);
#endif
            answers[i] += CN[p][q-h];
        }
    };

    auto addEdges = [&](std::vector<ui> & edges, ui a, ui b) {
        if(a > b) std::swap(a, b);

        auto st = pOEdge.begin() + pOIdx[a];
        auto ed = pOEdge.begin() + pOIdx[a + 1];
        ui idx = std::lower_bound(st, ed, b) - pOEdge.begin(); 

        edges.push_back(idx);
    };

    if(h == q) {
        for(auto i : HH) answers[i] += 1;
        return;
    }

    if(edC == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        updateAns(); return;
        return;
    }
    if(edC == 1) {
for(ui i = 0; i < p; i++) if(g.connectHash(C[0], Pn[i])) addEdges(PP, C[0], Pn[i]);
for(ui i = 0; i < h; i++) 
    if(g.connectHash(C[0], H[i])) 
        addEdges(PH, C[0], H[i]);
Pn[p++] = C[0];
    
        updateAns(); 

p -= 1;
for(ui i = 0; i < p; i++) if(g.connectHash(C[0], Pn[i])) PP.pop_back();
for(ui i = 0; i < h; i++)
    if(g.connectHash(C[0], H[i])) 
        PH.pop_back();
        return;
    }

    if(p + h + edC < q) return;

    // if(h == Q) {
    //     answers[Q] += 1;
    //     return;
    // }

    //init deg from the pre level
    if(sg.deg[deep].size() == 0) sg.deg[deep].resize(g.n);
    auto computeDegInC = [&](ui v) {
        ui costC = edC*2;
        ui costAdj = sg.degree2(v);
        ui costNAdj = nadj[v].size();
        ui d = 0;
        if(costC < costAdj && costC < costNAdj) {
            for(ui i = 0; i < edC; i++) if(g.connectHash(v, C[i])) d++;
        }
        else if(costAdj < costC && costAdj < costNAdj) {
            for(ui i = sg.pIdx[v]; i < sg.pIdx2[v]; i++) {
                ui w = sg.pEdge[i];
                if(C.isIn(w, 0, edC)) d++;
            }
        }
        else {
            d = edC;
            for(ui u : nadj[v]) if(C.isIn(u, 0, edC)) d--;
        }
        sg.deg[deep][v] = d;
    };
    for(ui i = 0; i < edC; i++) computeDegInC(C[i]);
    for(ui i = 0; i < edP; i++) computeDegInC(P[i]);

    //is P+C a k-plex, minDPC is the min degree in P+C
    ui minDPC = sg.deg[deep][C[0]] + edP - nn.getCntNonNei(C[0]);
    //find pivot in C with max degree in C
    ui maxD = sg.deg[deep][C[0]], pv = C[0];
    bool isPivot = false;
// ui nonNeiOfPv = nn.getCntNonNei(C[0]);
    for(ui i = 1; i < edC; i++) {
        ui v = C[i];
        ui d = sg.degreeDeep(deep, v);
        
        minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
        ui nonNei = nn.getCntNonNei(v);

        bool canServeAsPivot = true;
        for(ui j = 0; j < nonNei; j++) {
            ui w = nn.buffer[v*k + j];
            ui degWInC = sg.deg[deep][w];
            if(edC - degWInC + nn.getCntNonNei(w) > k) {
                canServeAsPivot = false;
                break;
            }
        }

        if(canServeAsPivot) {
#ifdef DDDEBUG
printf("canServeAsP : %u\n", v);
#endif
            if(isPivot) {
                if(d > maxD) {
                    maxD = d;
                    pv = v;
                }
            }
            else {
                isPivot = true;
                maxD = d;
                pv = v;
            }
        }
        else if(isPivot == false && d > maxD){
            maxD = d;
            pv = v;
            // nonNeiOfPv = nonNei;
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
        for(ui i = 0; i < edC; i++) {
            ui u = C[i];
            for(ui j = 0; j < p; j++) {
                if(g.connectHash(u, Pn[j]))
                    addEdges(PP, u, Pn[j]);
            }
            for(ui j = 0; j < h; j++) 
                if(g.connectHash(u, H[j])) 
                    addEdges(PH, u, H[j]);
            Pn[p++] = u;
        }
        
        bkPivot(deep + 1, 0, 0, p, h);

        for(ui i = edC-1; i >= 0; i--) {
            ui u = C[i];
            p--;
            for(ui j = 0; j < p; j++) if(g.connectHash(u, Pn[j])) PP.pop_back();
            for(ui j = 0; j < h; j++)
                if(g.connectHash(u, H[j])) 
                    PH.pop_back();
            if(i == 0) break;
        }
        return;
    }

    //build cand set by pv
    std::vector<ui> cand(edC);
    ui candSize = 0;
    
    ui countOfNonNeihborOfPv = nn.getCntNonNei(pv);
    ui C1Size = 0;
    for(ui j = 0; j < edC; j++) { 
        ui v = C[j];
        if(v == pv) continue;
        if(g.connectHash(pv, v)) {
            C.changeToByPos(j, C1Size++);
        }
        else {
            cand[candSize++] = v;
        }
    }
    if(!isPivot) cand[candSize++] = pv;
#ifdef DDDEBUG
printf("pv: %u", pv);
printf("cand:");
for(ui i = 0; i < candSize; i++) printf("%u ", cand[i]);printf("\n");
#endif

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

                newEnd = tmp;
            }
        }

        return newEnd;
    };

    auto solvePre = [&](ui v) {
        P.changeTo(v, edP++);
        C.changeTo(v, --edC);

        // for(ui j = 0; j < edP; j++) if(!g.connectHash(v, P[j])) nn.addNonNei(P[j], v);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.addNonNei(nn.buffer[v*k + i], v);
        }
        nn.addNonNei(v, v);

        //build tmpC for the next level
        //if v in P has k non-neighors
        ui tmpC = updateC(v, C, edC);

        if(tmpC*2 < nadj[v].size()) {
            for(ui j = 0; j < tmpC; ) {
                if(!g.connectHash(v, C[j])) {
                    if(nn.getCntNonNei(C[j])+1 >=k) {
                        C.changeToByPos(j, --tmpC);
                        continue;
                    }
                    nn.addNonNei(C[j], v);
                    j++;
                }
                else j++;
            }
        }
        else {
            for(ui u : nadj[v]) {
                if(C.isIn(u, 0, tmpC)) {
                    if(nn.getCntNonNei(u)+1 >= k) {
                        C.changeTo(u, --tmpC);
                    }
                    else nn.addNonNei(u, v);
                }
            }
        }

        return tmpC;
    };

    auto solveBack = [&](ui v, ui newEdC) {
        // for(ui j = 0; j < edP; j++) if(!g.connectHash(pv, P[j])) nn.pop(P[j]);
        nn.pop(v);
        for(ui i = 0, cntNN = nn.getCntNonNei(v); i < cntNN; i++) {
            nn.pop(nn.buffer[v*k + i]);
        }
        if(newEdC*2 < nadj[v].size()) {
            for(ui j = 0; j < newEdC; j++) 
                if(!g.connectHash(v, C[j])) nn.pop(C[j]);
        }
        else  for(ui u : nadj[v]) if(C.isIn(u, 0, newEdC)) nn.pop(u);
        // for(ui j = 0; j < edC; j++) if(!g.connectHash(v, C[j])) nn.pop(C[j]);
        P.changeTo(pv, --edP);
    };

    
//和pv没有公共非邻居的点，不用枚举
    
    if(isPivot) {
#ifdef DDDEBUG
printf("C2size == 0\n");
#endif

        P.changeTo(pv, edP++);
Pn[p] = pv;
for(ui i = 0; i < p; i++) if(g.connectHash(pv, Pn[i])) addEdges(PP, pv, Pn[i]);
for(ui i = 0; i < h; i++) if(g.connectHash(pv, H[i])) addEdges(PH, pv, H[i]);
        bkPivot(deep + 1, C1Size, edP, p+1, h);
        P.changeTo(pv, --edP);
for(ui i = 0; i < p; i++) if(g.connectHash(pv, Pn[i])) PP.pop_back();
for(ui i = 0; i < h; i++) if(g.connectHash(pv, H[i]))  PH.pop_back();
    }
    else{
#ifdef DDDEBUG
printf("C2size != 0\n");
#endif
        bkPivot(deep + 1, C1Size, edP, p, h);
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
        ui up = edP + k - nn.getCntNonNei(v) + bucket[0].size();
        for(ui j = 1; j < k; j++) {
            for(auto u: bucket[j]) {
                ui minJ = 0, minS = k + 1;
                // for(ui l = 0; l < edP; l++) {
                //     if(!g.connectHash(u, P[l])) {
                //         if(support[l] < minS) {
                //             minJ = l;
                //             minS = support[l];
                //         }
                //     }
                // }
                for(ui ll = 0, ct = nn.getCntNonNei(u); ll < ct; ll++) {
                    ui l = P.idx(nn.buffer[u*k+ll]);
                    if(support[l] < minS) {
                        minJ = l;
                        minS = support[l];
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
fflush(stdout);
#endif   
        if(nn.getCntNonNei(v) >= k) {
#ifdef DDDEBUG
printf("continue %u\n", v);
#endif
            C.changeTo(v, --edC);
            continue;
        }

        if(upPrune(v)) {
#ifdef DDDEBUG
printf("upPrune %u\n", v);
#endif
            C.changeTo(v, --edC);
            continue;
        }
H[h] = v;
for(ui j = 0; j < h; j++) if(g.connectHash(v, H[j]))  addEdges(HH, v, H[j]);
for(ui j = 0; j < p; j++) if(g.connectHash(v, Pn[j])) addEdges(PH, v, Pn[j]); 
        ui newEdC = solvePre(v);
        bkPivot(deep + 1, newEdC, edP, p, h + 1);
        solveBack(v, newEdC);
for(ui i = 0; i < h; i++) if(g.connectHash(v, H[i]))  HH.pop_back();
for(ui i = 0; i < p; i++) if(g.connectHash(v, Pn[i])) PH.pop_back();
    }
}


plexLocal::plexLocal(Graph && g, ui k, ui q):g(g), k(k), q(q) { 
    C.resize(g.n);
    P.resize(g.n);

    sg.pIdx.resize(g.n);
    sg.pIdx2.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(g.n);
    sg.deg[0].resize(g.n);

    nn.init(g.n, k);

    bucket.resize(k);
    support.resize(g.n);

    answers.resize(g.m);

    maxSize = g.coreNumber*g.coreNumber;
    maxD = maxSize;
    maxD2 = 20;
    computeC();

    pOEdge.resize(g.m / 2);
    pOIdx.resize(g.n + 1);
    for(ui u = 0; u < g.n; u++) {
        pOIdx[u + 1] = pOIdx[u];
        for(ui j = g.pIdx2[u]; j < g.pIdx[u + 1]; j++) {
            pOEdge[pOIdx[u+1]++] = g.pEdge[j];
        }
    }

    // P, H, PP, PH, HH
    Pn.resize(maxSize);
    H.resize(maxSize);

    ui maxSize2 = std::min(g.m/2,maxSize * maxSize);
    PP.resize(maxSize2);
    PH.resize(maxSize2);
    HH.resize(maxSize2);
    PP.clear();
    PH.clear();
    HH.clear();

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

    printf("plexELocal.h\n");
}