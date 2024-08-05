#include "kPList.h"
#include <algorithm>

// #define BASELINE
// #define DDDEBUG

// #define PRINT

#ifdef DDDEBUG
#include <iostream>
#endif

ull kplist::run() {
#ifdef DDDEBUG
g.print();
#endif
    printf("kPList v2::run\n");

    g.initHash();
    printf("init Hash\n");

    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<' '<<answer<<std::endl;
#endif
// if(u % 1000 == 0) {
//     printf("%u\n", u);fflush(stdout);
// }
        //P is empty
        //X is 2-hop neighbors < u
        //C is 2-hop neighbors > u
        if(g.pIdx[u + 1] - g.pIdx2[u] + k < q) continue;
        std::vector<ui> C;
        // P.clear(); 
        edP = 0;
        C.clear();

        // P.push_back(u);
        P.changeTo(u, edP++);
        level[u] = 1;
        nn.addNonNei(u, u);

        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.push_back(g.pEdge[i]);
            level[g.pEdge[i]] = 1;
        }
#ifdef DDDEBUG
printf("C:"); for(ui i = 0; i < C.size(); i++) printf("%u ", C[i]); printf("\n");
printf("edC %u\n", C.size());
#endif
        //delete C
        while(true) {
            ui sz = C.size();
            std::vector<ui> C2;
            std::vector<ui> C3;
            for(ui i = 0; i < sz; i++) {
                ui v = C[i];
                ui d = 0;
                if(g.degree(v) <= sz) {
                    for(ui i = g.pIdx[v]; i < g.pIdx[v + 1]; i++) {
                        if(level[g.pEdge[i]] == 1) {
                            d++;
                            if(d + 2*k >= q) break;
                        }
                    }
                }
                else {
                    for(ui j = 0; j < sz; j++) {
                        ui w = C[j];

                        if(g.connectHash(v, w)) {
                            d++;
                            if(d + 2*k >= q) break;
                        }
                    }
                }

                if(d + 2*k >= q) C2.push_back(v);
                else C3.push_back(v);
            }

            if(C2.size() == sz) break;
            C = C2;
            for(auto v:C3) level[v] = 0;
        }
#ifdef DDDEBUG
printf("After delete C\n");
printf("C:"); for(ui i = 0; i < C.size(); i++) printf("%u ", C[i]); printf("\n");
#endif
        //two hop neighbors
        for(ui j = 0, ed = C.size(); j < ed; j++) {
            ui v = C[j];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(w <= u) continue;
                if(g.connectHash(u, w)) continue;
                if(level[w] == 1) continue; 
                ui d = 0;
                for(ui l = 0; l < ed; l++) {
                    ui x = C[l]; 
                    if(g.connectHash(w, x)) {
                        d++;
                        if(d + 2 * k >= q + 2) {
#ifdef DDDEBUG
printf("Add 2-hop neighbor %u of %u in C\n", w, v);
#endif
                            C.push_back(w);
// assert(nn.getCntNonNei(w) == 0);
                            nn.addNonNei(w, u);
                            level[w] = 1;
                            break;
                        }
                    }
                }
            }
        }

        std::sort(C.begin(), C.end());
        listing(1, C);

        nn.clear(u);
        for(ui i = 0; i < C.size(); i++) nn.clear(C[i]);
        level[u] = 0;
        for(ui i = 0; i < C.size(); i++) level[C[i]] = 0;
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
    return answer;
}

void kplist::listing(ui deep, const std::vector<ui> & C) {
#ifdef DDDEBUG
printf("    DDDDDeep %u\n", deep);
printf("P:");
for(auto v:P) printf("%u ", v);printf("\n");
printf("C:");
for(auto v:C) printf("%u ", v);printf("\n");
printf("NonNeis:\n");
for(auto v:P) {
    printf("%u:", v);
    for(ui i = 0; i < nn.getCntNonNei(v); i++) 
        printf("%u ", nn.buffer[v*k+i]);
    printf("\n");
}
for(auto v:C) {
    printf("%u:", v);
    for(ui i = 0; i < nn.getCntNonNei(v); i++) 
        printf("%u ", nn.buffer[v*k+i]);
    printf("\n");
}
#endif

// if(deep == 1) {
//     for(ui i = 1; i < C.size(); i++) {
//         assert(C[i] > C[i-1]);
//     }
// }
// if(P.size() == k) {
// if(C.size() > k*g.coreNumber) {
// printf("%u %u\n", C.size(), g.coreNumber*k);fflush(stdout);
// for(ui u : P) {
//     ui d = 0;
//     for(auto v: C) if(g.connectHash(u, v)) {
//         assert(u < v);
//         d++;
//     }
//     assert(d <= g.pIdx[u + 1] - g.pIdx2[u]);
// }
// // for(ui i = 0; i < C.size(); i++) assert(nn.getCntNonNei(C[i]) < k);
// }
//     assert(C.size() <= k*g.coreNumber);
// }

    if(edP == q - 1) {
#ifdef DDDEBUG
printf("ans+%u\n", C.size());
#endif

#ifdef PRINT
for(auto v: C) {
    for(auto u:P) printf("%u ", u);
    printf(" %u\n", v);
}
#endif
        answer += C.size();
        return;
    }

    if(C.size() + edP < q) {
#ifdef DDDEBUG
printf("C+P<q, %u %u %u\n", C.size(), P.size(), q);
#endif
        return; 
    }

    auto upPrune = [&](ui i, ui v) {
        //the prune technique of Dai
        for(ui j = 0; j < k; j++) bucket[j].clear();
     
        for(ui j = i + 1; j < C.size(); j++) {
            if(g.connectHash(v, C[j])) {
                ui cNei = nn.getCntNonNei(C[j]);
                bucket[cNei].push_back(C[j]);  
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
                // for(ui l = 0; l < P.size(); l++) {
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

    for(ui i = 0; i < C.size(); i++) {
        std::vector<ui> newC;
        ui u = C[i];

        if(upPrune(i, u)) continue;

        // P.push_back(u);
        P.changeTo(u, edP++);
        level[u] = deep + 1;
        nn.addNonNei(u, u);
#ifdef DDDEBUG
printf("    deep %u, cand %u\n", deep, u);
#endif

        for(ui l = i + 1; l < C.size(); l++) level[C[l]] = deep + 1;
        for(ui j = 0, cntNN = nn.getCntNonNei(u); j < cntNN; j++) {
            ui w = nn.buffer[u*k + j];
            if(w != u) nn.addNonNei(w, u);

            if(nn.getCntNonNei(w) < k) continue;

            for(ui l = i + 1; l < C.size(); l++) {
                if(!g.connectHash(w, C[l])) {
                    level[C[l]] = deep;
                }
            }
        }
        
        for(ui j = i + 1; j < C.size(); j++) {
            if(level[C[j]] != deep + 1) continue;

            if(g.connectHash(u, C[j])) {
                newC.push_back(C[j]);
            }
            else {
                if(nn.getCntNonNei(C[j]) + 1 >= k) {
                    level[C[j]] = deep;
                    continue;
                }
                nn.addNonNei(C[j], u);
                newC.push_back(C[j]);
            }
        }

// auto isKPlex = [&](ui v) -> bool {
//     ui nd = 1;
//     for(auto w:P) {
//         if(!g.connectHash(v, w)) {
//             nd++;
//             if(nn.getCntNonNei(w) == k) {
// printf("tp1, w%u u%u v%u wu%d, vw%d wv%d %u\n", 
//     w, u, v, (int)g.connectHash(w, u), (int)g.connectHash(v, w),
//     g.connectHash(w, v), deep);fflush(stdout);
//                 return false;
//             }
//         }
//     }
//     if(nd > k) return false;
//     return true;
// };
// for(auto v : newC) assert(isKPlex(v));
        
        if(newC.size()) listing(deep + 1, newC);

        // P.pop_back();
        edP--;
        // nn.pop(u);
        level[u] = deep;
        for(ui j = 0, cntNN = nn.getCntNonNei(u); j < cntNN; j++) {
            ui w = nn.buffer[u*k + j];
            nn.pop(w);
        }
        for(auto v : newC) if(!g.connectHash(u, v)) nn.pop(v);
        for(auto v : newC) level[v] = deep;
    }
}