#include "plexCounting.h"
#include <iostream>

// #define BASELINE

// #define DDDEBUG

// #define C1C2
#ifdef C1C2
ui c1 = 0, c2 = 0;
#endif

std::vector<ull> plexCounting::run() {
#ifdef DDDEBUG
g.print();
#endif
    printf("plexCountingV5.cpp::run (not right algirthm)\n");
    
    std::vector<ui> toDelete(std::min(g.n, g.maxD * k));
    toDelete.clear();

    g.initHash();
    printf("init Hash\n");

    // Pnon.resize(std::min(g.n, g.maxD * k + 1));
    // Pnon.clear();
    // H.resize(std::min(g.n, g.maxD * k + 1));
    // H.clear();
    // Pr.resize(g.n);
    // dp.resize(std::min(g.n, g.maxD * k + 1));
    Pnon.resize(g.coreNumber*g.coreNumber);
    Pnon.clear();
    H.resize(g.coreNumber*g.coreNumber);
    H.clear();
    Pr.resize(g.n);
    Pr.clear();
    dp.resize(g.coreNumber*g.coreNumber);
    dp[0] = 1.0;

    nadj.resize(g.n);

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
#ifdef DDDEBUG
double pa = answers[q];
#endif
        H.push_back(u);
        bkPivot(1, edC, edP, 0, 1);
        H.clear();
#ifdef DDDEBUG
printf("diff %u: %.0f\n", u, answers[q]-pa);
#endif
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
    }
}
printf("cnt:%u\n", cnt);

#endif

#ifdef C1C2
printf("c1:%u\n", c1);
printf("c2:%u\n", c2);
#endif
    return answers;
}


void plexCounting::bkPivot(ui deep, ui edC, ui edP, ui p, ui h) {
#ifdef DDDEBUG
printf("       deep %u, p %u, h %u\n", deep, p, h);
printf("C:"); for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:"); for(ui i = 0; i < edP; i++) printf("%u ", P[i]); printf("\n");
printf("Pnon:"); for(ui i = 0; i < Pnon.size(); i++) printf("%u ", Pnon[i]); printf("\n");
printf("H:"); for(ui i = 0; i < H.size(); i++) printf("%u ", H[i]); printf("\n");
printf("\n");
#endif


    if(edC == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        for(auto w : H) Pr[w] = 0;
        
        for(auto v : Pnon) {
            ui w = nn.buffer[v*k];
            if(nn.getCntNonNei(w) < k) Pr[w]++;
        }
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        ui sum = 0;
        for(auto w : H) if(nn.getCntNonNei(w) < k) 
            sum += std::min(Pr[w], k-nn.getCntNonNei(w));

        for(ui i = 1; i < sum; i++) dp[i] = 0.0;

        // dp[0] = 1.0;
        ui presum = 0;
        for(auto w : H) if(nn.getCntNonNei(w) < k) {
            ui t = std::min(Pr[w], k-nn.getCntNonNei(w));
            for(ui j = presum; j >= 0; j--) {
                for(ui x = 1; x <= t; x++) {
                    dp[j+x] += dp[j]*CN[Pr[w]][x];
                }
                if(j == 0) break;
            }
            presum += t;
        }

        for(ui j = 0; j <= presum; j++) {
            if(dp[j] < 0.5) break;
            for(ui i = 0; i <= p; i++) if(h+i+j <= Q) {
#ifdef DDDEBUG
if(h+i+j==q)
printf("ans+%.0f\n", dp[j] * CN[p][i]);
#endif
                answers[h+i+j] += dp[j] * CN[p][i];
            }
        }

        return;
    }

    if(p + h + edC + Pnon.size() < q) return;
    if(h == Q-1) {
        answers[Q-1] += 1;
        ui pr = 0;
        for(auto v : Pnon) {
            ui w = nn.buffer[v*k];
            if(nn.getCntNonNei(w) < k) pr++;
        }
        answers[Q] += p + edC + pr;
        return;
    }
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
    // for(ui i = 0; i < edP; i++) computeDegInC(P[i]);
// #ifdef DDDEBUG
// printf("CdeginC:");
// for(ui i = 0; i < edC; i++) printf("%u-%u ", C[i], sg.deg[deep][C[i]]);
// printf("\nPdeginC:");
// for(ui i = 0; i < edP; i++) printf("%u-%u ", P[i], sg.deg[deep][P[i]]);
// printf("\n");
// #endif

    //is P+C a k-plex, minDPC is the min degree in P+C
    // ui minDPC = sg.deg[deep][C[0]] + edP - nn.getCntNonNei(C[0]);
    //find pivot in C with max degree in C
    ui maxD = sg.deg[deep][C[0]], pv = C[0];
ui nonNeiOfPv = nn.getCntNonNei(C[0]);
    for(ui i = 1; i < edC; i++) {
        ui v = C[i];
        ui d = sg.degreeDeep(deep, v);
        
        // minDPC = std::min(minDPC, d + edP - nn.getCntNonNei(v));
ui nonNei = nn.getCntNonNei(v);

if(nonNeiOfPv <= 1) {
    if(d > maxD) {
        maxD = d;
        pv = v;
        nonNeiOfPv = nonNei;
    }
}
    }

    if(nonNeiOfPv > 1) {
        maxD = sg.deg[deep][C[0]];
        pv = C[0];
        for(ui i = 1; i < edC; i++) {
            ui v = C[i];
            ui d = sg.degreeDeep(deep, v);
            if(d > maxD) {
                maxD = d;
                pv = v;
            }
        }
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
            // bool isInC1 = true;
            for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
                ui w = nn.buffer[pv*k + i];
                if(!g.connectHash(v, w)) {
                    // cand[candSize++] = v;
                    C2Size++;
                    // isInC1 = false;
                    break;
                }
            }
            // if(isInC1) 
            C.changeToByPos(j, C1Size++);
            // if(isInC1) C1Size++;
        }
        else {
            // if(canAdd(v)) C2Size++;
            // if(countOfNonNeihborOfPv+1 < k && nn.getCntNonNei(v)+1 < k) {
            //     bool canAdd = true;
            //     for(ui i = 0; i < countOfNonNeihborOfPv; i++) {
            //         ui w = nn.buffer[pv*k + i];
            //         if(nn.getCntNonNei(w)+1 == k && !g.connectHash(v, w)) {
            //             canAdd = false;
            //         }
            //     }
            //     if(canAdd) C2Size++;
            // }
            
            cand[candSize++] = v;
        }
    }
    if(C2Size > 0 && nonNeiOfPv > 1) cand[candSize++] = pv;
#ifdef DDDEBUG
printf("pv: %u\n", pv);
printf("cand:");
for(ui i = 0; i < candSize; i++) printf("%u ", cand[i]);
printf("\n");
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

    

#ifdef C1C2
c2++;
#endif 
//和pv没有公共非邻居的点，不用枚举
    
    if(C2Size == 0 && nonNeiOfPv == 0) {
#ifdef DDDEBUG
printf("C2size==0, deep %u, pv %u\n", deep, pv);
#endif 
        P.changeTo(pv, edP++);
        bkPivot(deep + 1, C1Size, edP, p+1, h);
        P.changeTo(pv, --edP);
    }
    else if(nonNeiOfPv == 1) {//if nonNeiOfPv==0, C2Size==0
#ifdef DDDEBUG
printf("nonNeiOfPv==1, deep %u, pv %u\n", deep, pv);
#endif 
        Pnon.push_back(pv);
        bkPivot(deep + 1, C1Size, edP, p, h);
        Pnon.pop_back();
    }
    else {
#ifdef DDDEBUG
printf("listing, deep %u, pv %u\n", deep, pv);
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
            continue;
        }

        if(upPrune(v)) {
#ifdef DDDEBUG
printf("upPrune %u\n", v);
#endif
            C.changeTo(v, --edC);
            continue;
        }

        ui newEdC = solvePre(v);
        H.push_back(v);
        bkPivot(deep + 1, newEdC, edP, p, h + 1);
        H.pop_back();
        solveBack(v, newEdC);
    }
}
