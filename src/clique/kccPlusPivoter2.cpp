#include "kccPlusPivoter.h"
#include <cassert>

// #define DDEBUG
// #define DEBUG
// #define BASELINE
// #define BRANCHES

#ifdef BRANCHES
ui branches = 0;
#endif

#ifdef DDEBUG
ui uu = 0, uuu = 2;
#endif
double kccPlusPivoter::run() {
    printf("kccPlusPivoter2.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
if(uu==uuu)
printf("        st %u\n", u);
#endif
// if(u==188) {
// // if(u==1696186) {//ski
// for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
//     ui v = g.pEdge[i];
//     for(ui j = i + 1; j < g.pIdx[u + 1]; j++) {
//         ui w = g.pEdge[j];
//         if(g.connectHash(v, w)) printf("%u %u\n", v, w);
//     }
// }
// }
// else continue;
        ui deep = 0;
        nodes[deep].clear();
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            nodes[deep].push_back(g.pEdge[i]);
            // level[g.pEdge[i]] = deep;
        }

#ifdef DDEBUG
double p = answer;
#endif
        if(nodes[deep].size() == 0) {
// printf("ut %u 0\n", u);
            continue;
        }
                
        std::vector<ui> & C = nodes[deep];
        std::vector<ui> & nxtC = nodes[deep + 1];
        nxtC.clear();
        ui maxDeg = 0, maxV = C[0];
        auto findM = [&](std::vector<ui> & C, ui & maxDeg, ui & maxV) {
            for(ui i = 0; i < C.size(); i++) {
                ui deg = 0;
                for(ui j = 0; j < C.size(); j++) {
                    if(g.connectHash(C[i], C[j])) deg++;
                }
                if(deg > maxDeg) {
                    maxDeg = deg;
                    maxV = C[i];
                }
            }
        };
        findM(C, maxDeg, maxV);
        std::vector<ui> cliqueNei;
        clique.clear(); cliqueNei.clear();
        clique.push_back(maxV);
        for(auto v : C) {
            if(g.connectHash(maxV, v)) cliqueNei.push_back(v);
            else if(v != maxV) nxtC.push_back(v);
        }
// assert(C.size() == g.pIdx[u+1] - g.pIdx2[u]);
// if(nxtC.size() + clique.size() + cliqueNei.size() != C.size()) {
//     printf("%u %u %u %u\n", nxtC.size(), clique.size(), cliqueNei.size(), C.size());
//     fflush(stdout);
// }
// assert(nxtC.size() + clique.size() + cliqueNei.size() == C.size());

        while(cliqueNei.size() > 0) {
            ui maxDeg = 0, maxV = cliqueNei[0];
            findM(cliqueNei, maxDeg, maxV);
            clique.push_back(maxV);

            std::vector<ui> newCliqueNei;
            for(ui i = 0; i < cliqueNei.size(); i++) {
                if(g.connectHash(maxV, cliqueNei[i])) {
                    newCliqueNei.push_back(cliqueNei[i]);
                }
                else if(cliqueNei[i] != maxV) nxtC.push_back(cliqueNei[i]);
            }
            cliqueNei = newCliqueNei;
        }
// for(ui i = 0; i < clique.size(); i++) {
//     for(ui j = i + 1; j < clique.size(); j++) {
//         assert(g.connectHash(clique[i], clique[j]));
//     }
// }
// assert(nxtC.size() + clique.size() == C.size());
        deep += 1;
        for(ui v : nodes[deep]) sg.pIdx[v] = sg.deg[deep][v] = g.pIdx[v];
        for(ui v : nodes[deep]) level[v] = deep;
        for(ui v : nodes[deep]) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == deep) {
                    sg.pEdge[sg.deg[deep][v]++] = w;
                    sg.pEdge[sg.deg[deep][w]++] = v;
                }
            }
        }

        for(ui v : clique) degR[deep][v] = 0;
        for(ui v : nodes[deep]) degL[deep][v] = 0;
        for(ui v : clique) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == deep) {
                    R[v][degR[deep][v]++] = w;
                    L[w][degL[deep][w]++] = v;
                }
            }
        }

        for(ui v : nodes[deep]) level[v] = 0;
        for(ui v : clique) level[v] = deep;
        for(ui v : nodes[deep]) {
            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == deep) {
                    L[v][degL[deep][v]++] = w;
                    R[w][degR[deep][w]++] = v;
                }
            }
        }
        for(ui v : nodes[deep]) level[v] = deep;

#ifdef DDEBUG
if(u==uuu) {
    printf("subGraph:\n");
    for(ui i = 0; i < nodes[deep].size(); i++) {
        ui v = nodes[deep][i];
        printf("%u:", v);
        for(ui j = i + 1; j < nodes[deep].size(); j++) {
            ui w = nodes[deep][j];
            if(g.connectHash(v, w)) printf(" %u", w);
        }
        printf("\n");
    }
}
#endif
#ifdef DDEBUG
if(u==uuu) {
    printf("subGraphClique:\n");
    for(ui i = 0; i < clique.size(); i++) {
        ui v = clique[i];
        printf("%u:", v);
        for(ui j = 0; j < degR[deep][v]; j++) {
            ui w = R[v][j];
            printf(" %u", w);
        }
        printf("\n");
    }
}
#endif
        Ps.clear(); Hs.clear();
        Hs.push_back(u);

        searchSGClique(deep, 0, 1, clique.size());

// for(ui v : nodes[deep]) level[v] = deep;
        for(ui v : nodes[deep]) {
            ui & ed = sg.deg[deep][v];
            ui tmp = sg.deg[deep][v];
            ed = sg.pIdx[v];
            for(ui j = sg.pIdx[v]; j < tmp; j++) {
                if(sg.pEdge[j] > v) sg.pEdge[ed++] = sg.pEdge[j];
            }
        }
        listClique(deep, clique.size());

        for(ui v : nodes[deep]) level[v] = 0;
        for(ui v : clique) level[v] = 0;
// printf("ut %u %.0f\n", u, answers[k] - p);
#ifdef DDEBUG
printf("diff %u %.0f\n", u, answer - p);
#endif
    }

#ifdef BASELINE
std::vector<ui> vc(g.n, 0);
auto print = [&](uint32_t x) {
    printf("c:");
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        vc[u]++;
        break;
    }
};
auto check = [&](ui x) {
    ui sz = 0;
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) sz++;
    if(sz != k) return false;

    for(ui u = 0; u < g.n; u++) if((1<<u) & x){
        ui d = 0;
        for(ui v = 0; v < g.n; v++) if((1<<v) & x){
            if(g.connectHash(u, v)) d++;
        }
        if(d+1 < sz) return false; 
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
for(ui u = 0; u < g.n; u++) printf("u_cntu %u %u\n", u, vc[u]);
#endif
#ifdef BRANCHES
printf("branches:%u\n", branches);
#endif
    return answer;
}

void kccPlusPivoter::searchSGClique(ui deep, ui p, ui h, ui edClique) {
#ifdef BRANCHES
branches++;
#endif
#ifdef DDEBUG
if(uu == uuu) {
printf("    deep %u\nC:", deep);
for(ui i = 0; i < nodes[deep].size(); i++) {
printf("%u ", nodes[deep][i]);
}printf("\n");
printf("p %u, h %u\nclique%u:", p, h, edClique);
for(ui i = 0; i < edClique; i++) {
printf("%u ", clique[i]);
}printf("\n");
printf("R:");
for(ui i = 0; i < edClique; i++) {
    printf("%u:", clique[i]);
    for(ui j = 0; j < degR[deep][clique[i]]; j++) printf("%u ", R[clique[i]][j]);
}printf("\n");
printf("L:");
for(ui i = 0; i < Ps.size(); i++) {
    printf("%u:", Ps[i]);
    for(ui j = 0; j < degL[deep][Ps[i]]; j++) printf("%u ", L[Ps[i]][j]);
}
for(ui i = 0; i < nodes[deep].size(); i++) {
    printf("%u:", nodes[deep][i]);
    for(ui j = 0; j < degL[deep][nodes[deep][i]]; j++) printf("%u ", L[nodes[deep][i]][j]);
}printf("\n");
}
#endif

    if(h == k-1) {
        answer += p + nodes[deep].size() + edClique;
#ifdef DDEBUG
if(uu == uuu) {
    printf("h==k-1, ans+=%u+%u+%u\n", p, nodes[deep].size(), edClique);
}
#endif
        return;
    }
    if(p+h+nodes[deep].size()+2 < k) {
#ifdef DDEBUG
if(uu == uuu) {
    printf("return 2\n");
}
#endif
        return;
    }
    if(nodes[deep].size() + p + h + edClique < k) {
#ifdef DDEBUG
if(uu == uuu) {
    printf("return 3\n");
}
#endif
        return;
    }
    
    //choose 0,1,2 vertices from clique
    auto updateAns = [&]() {
        // return;
#ifdef DDEBUG
if(uu == uuu) {
printf("\n In updateAns\n");
printf("deep %u\nC:", deep);
for(ui i = 0; i < nodes[deep].size(); i++) {
printf("%u ", nodes[deep][i]);
}printf("\n");
printf("p %u, h %u\nclique:", p, h);
for(ui i = 0; i < edClique; i++) {
printf("%u ", clique[i]);
}printf("\n");
printf("Ps:");
for(ui i = 0; i < p; i++) printf("%u ", Ps[i]); printf("\n");
printf("Hs:");
for(ui i = 0; i < Hs.size(); i++) printf("%u ", Hs[i]); printf("\n");
printf("R:");
for(ui i = 0; i < edClique; i++) {
    printf("%u:", clique[i]);
    for(ui j = 0; j < degR[deep][clique[i]]; j++) printf("%u ", R[clique[i]][j]);
}printf("\n");
printf("L:");
for(ui i = 0; i < Ps.size(); i++) {
    printf("%u:", Ps[i]);
    for(ui j = 0; j < degL[deep][Ps[i]]; j++) printf("%u ", L[Ps[i]][j]);
}printf("\n");
}
#endif
        
        //choose 0 vertices from clique
        if(p >= k - h) answer += CN[p][k - h];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans0 += C(%u,%u)\n", p, k-h);
}
#endif
        //choose 1 vertices from clique
        for(ui i = 0; i < edClique; i++) {
            if(degR[deep][clique[i]] >= k-h-1) {
                answer += CN[degR[deep][clique[i]]][k-h-1];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans1 += C(%u,%u), clique node %u\n", degR[deep][clique[i]] , k-h-1, clique[i]);
}
#endif
            }
        }

        //choose 2 vertices from clique
        if(k < h+2) return;
        if(k == h+2) {
            answer += CN[edClique][2];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans1.5 += C(%u,%u)\n", edClique, 2);
}
#endif
            return;
        }
        for(ui i = 0; i < edClique; i++) weight[clique[i]] = 0;
        twoHopNodes.clear();
        for(ui i = 0; i < edClique; i++) {
            ui u = clique[i];
            for(ui j = 0; j < degR[deep][u]; j++) {
                ui v = R[u][j];

                for(int l = 0; l < degL[deep][v]; l++) {
                    ui w = L[v][l];
#ifdef DDEBUG
if(uu == uuu) {
    printf("%u-%u-%u\n", u, v, w);
}
#endif      
                    if(w > u) {
                        weight[w]++;
                        if(weight[w] == 1) {
                            twoHopNodes.push_back(w);
#ifdef DDEBUG
if(uu == uuu) {
    printf("pu 2 hop %u\n", w);
}
#endif    
                        }
                    }
                }
            }

            for(auto v : twoHopNodes) {
                if(weight[v] >= k-h-2) {
                    answer += CN[weight[v]][k-h-2];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans2 += C(%u,%u), clique node %u-%u\n", weight[v], k-h-2, clique[i], v);
}
#endif
                }
                weight[v] = 0;
            }
        }
        
    };
    std::vector<ui> & C = nodes[deep];
    ui edC = C.size();

    if(edC == 0) {
#ifdef DDEBUG
if(uu == uuu) {
printf("edC==0\n");
}
#endif
        updateAns(); return;
    }
    if(edC == 1) {
#ifdef DDEBUG
if(uu == uuu) {
printf("edC==1\n");
}
#endif
        Ps.push_back(C[0]);
        p += 1; updateAns(); 
        Ps.pop_back();
        return;
    }
    if(edC == 2) {
        if(sg.deg[deep][C[0]] > sg.pIdx[C[0]]) {
#ifdef DDEBUG
if(uu == uuu) {
printf("edC==2 and connect\n");
}
#endif
            Ps.push_back(C[0]);
            Ps.push_back(C[1]);
            p += 2; updateAns();
            Ps.pop_back();
            Ps.pop_back();
        }
        else {
#ifdef DDEBUG
if(uu == uuu) {
printf("edC==2 not connected, 1\n");
}
#endif
            Ps.push_back(C[0]);
            p += 1;
            for(ui i = 0; i < degL[deep][C[1]]; i++) {
                ui v = L[C[1]][i];
                for(ui j = 0; j < degR[deep][v]; j++) {
                    if(R[v][j] == C[1]) {
                        std::swap(R[v][j], R[v][--degR[deep][v]]);
                        break;
                    }
                }
            }
            updateAns();
            // for(ui i = 0; i < degL[deep][C[1]]; i++) {
            //     ui v = L[C[1]][i];
            //     degR[deep][v]++;
            // }

            Ps.pop_back();
            Hs.push_back(C[1]);
            p -= 1; h += 1; 
            for(ui j = 0; j < edClique; ) {
                ui u = clique[j];
                if(!g.connectHash(C[1], u)) {
                    level[u] = deep-1;
                    for(ui i = 0; i < degR[deep][u]; i++) {
                        ui v = R[u][i];
                        for(ui l = 0; l < degL[deep][v]; l++) {
                            if(L[v][l] == u) {
                                std::swap(L[v][l], L[v][--degL[deep][v]]);
                                break;
                            }
                        }
                    }
                    std::swap(clique[j], clique[--edClique]);
                }
                else {
// assert(level[clique[j]] == deep);
                    j++;
                }
            }

            for(ui i = 0; i < degL[deep][C[0]]; i++) {
                ui v = L[C[0]][i];
#ifdef DDEBUG
if(uu == uuu) {
printf("del %u-%u\n", C[0], v);
}
#endif
                if(level[v] != deep) continue;
                for(ui j = 0; j < degR[deep][v]; j++) {
                    if(R[v][j] == C[0]) {
                        std::swap(R[v][j], R[v][--degR[deep][v]]);
                        break;
                    }
                }
            }
// for(ui i = 0; i < edClique; i++) assert(degR[deep][clique[i]] == 0);
#ifdef DDEBUG
if(uu == uuu) {
printf("edC==2 not connected, 2\n");
}
#endif
            updateAns();
            Hs.pop_back();
        }
        return;
    }

    ui pivot = C[0], pivotDeg = 0, num = 0;
    for(ui i = 0; i < edC; i++) {
        // ui tmp = 0;
        ui v = C[i];
        ui tmp = sg.deg[deep][v] - sg.pIdx[v];

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }
#ifdef DDEBUG
if(uu == uuu) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    
    if(pivotDeg+1 == edC && num == edC) {
        p += edC; 
        for(ui i = 0; i < edC; i++) {
            Ps.push_back(C[i]);
            // level[C[i]] = deep+1;
        }
#ifdef DDEBUG
if(uu == uuu) {
printf("cand is a clique\n");
}
#endif
        updateAns();
        for(ui i = 0; i < edC; i++) {
            Ps.pop_back();
            // level[C[i]] = deep;
        }
        return;
    }

    std::vector<ui> & nxtC = nodes[deep + 1];
    nxtC.clear();
    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui w = sg.pEdge[i];
        nxtC.push_back(w);
        level[w] = deep + 1;
    }
    auto updateSG = [&]() {
        for(ui j = 0; j < nxtC.size(); j++) {
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
    updateSG();

    auto updateR = [&](ui edClique) {
        for(ui i = 0; i < edClique; i++) {
            ui u = clique[i];
            ui & ed = degR[deep + 1][u];
            ed = degR[deep][u];
            for(ui l = 0; l < ed; ) {
                ui w = R[u][l];
                if(level[w] != deep + 1 && !inP[w]) {
                    std::swap(R[u][l], R[u][--ed]);
                }
                else l++;
            }
        }
    };
    Ps.push_back(pivot);
    inP[pivot] = true;
    updateR(edClique);
    for(ui i = 0; i < edClique; i++) level[clique[i]] = deep + 1;
    for(auto v : nxtC) degL[deep+1][v] = degL[deep][v];
    for(auto v : Ps) degL[deep+1][v] = degL[deep][v];

    searchSGClique(deep + 1, p + 1, h, edClique);
    
    Ps.pop_back();
    inP[pivot] = false;
    for(auto v:nxtC) level[v] = deep;
    for(ui i = 0; i < edClique; i++) level[clique[i]] = deep;

#ifdef DDEBUG
if(uu == uuu) {
printf("    redeep %u, \nC:", deep);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("R:");
for(ui i = 0; i < edClique; i++) {
    printf("%u:", clique[i]);
    for(ui j = 0; j < degR[deep][clique[i]]; j++) printf("%u ", R[clique[i]][j]);
}printf("\n");
printf("L:");
for(ui i = 0; i < Ps.size(); i++) {
    printf("%u:", Ps[i]);
    for(ui j = 0; j < degL[deep][Ps[i]]; j++) printf("%u ", L[Ps[i]][j]);
}printf("\n");
}
for(ui i = 0; i < nodes[deep].size(); i++) {
    if(nodes[deep][i] == pivot) continue;
    // if(g.connectHash(pivot, nodes[deep][i])) continue;
    printf("%u:", nodes[deep][i]);
    for(ui j = 0; j < degL[deep][nodes[deep][i]]; j++) printf("%u ", L[nodes[deep][i]][j]);
}printf("\n");
#endif
    if(pivotDeg+1 < C.size())
    for(ui u : C) {
        if(u == pivot) continue;
        if(g.connectHash(pivot, u)) continue;
        level[u] = deep - 1;
#ifdef DDEBUG
if(uu == uuu) {
printf("    Cdeep %u, cand %u\nC:", deep, u);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
        nxtC.clear();
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui w = sg.pEdge[j];
            if(level[w] != deep) continue;
// assert(level[w] == deep);
            level[w] = deep + 1;
            nxtC.push_back(w);
        }
        updateSG();

        Hs.push_back(u);
        ui newEdClique = edClique;
        for(ui j = 0; j < newEdClique; ) {
            if(!g.connectHash(u, clique[j])) {
                std::swap(clique[j], clique[--newEdClique]);
            }
            else j++;
        }
        for(ui j = 0; j < newEdClique; j++) level[clique[j]] = deep+1;
        
        updateR(newEdClique);
        auto updateL = [&](ui v) {
            ui & ed = degL[deep + 1][v];
            ed = degL[deep][v];
            for(ui l = 0; l < ed; ) {
                ui w = L[v][l];
                if(level[w] != deep + 1) {
                    std::swap(L[v][l], L[v][--ed]);
                }
                else l++;
            }
        };
        for(ui v : nxtC) updateL(v);
        for(ui v : Ps) updateL(v);

        searchSGClique(deep + 1, p, h + 1, newEdClique);
        Hs.pop_back();

        for(auto v : nxtC) level[v] = deep;
        for(ui j = 0; j < newEdClique; j++) level[clique[j]] = deep;
    }
}

void kccPlusPivoter::listClique(ui deep, ui edClique) {
#ifdef DDEBUG
if(uu==uuu){
printf("    listingdeep %u\n", deep);
printf("C:");
for(auto v: nodes[deep]) printf("%u ", v);printf("\n");
printf("clique:");
for(ui i = 0; i < edClique; i++) printf("%u ", clique[i]);printf("\n");
}
#endif

    std::vector<ui> & C = nodes[deep];

    if(k - deep <= 3) {
        if(k - deep == 3 && edClique >= 3) {
            answer += CN[edClique][3];
        }
        return;
    }

    if(edClique >= k - deep) {
#ifdef DDEBUG
if(uu==uuu){
printf("ans += C(%u, %u)\n", edClique, k-deep);
}
#endif
        answer += CN[edClique][k - deep];
    }

    std::vector<ui> & nxtC = nodes[deep + 1];
    
    for(ui i = 0; i < C.size(); i++) {
        ui u = C[i];
#ifdef DDEBUG
if(uu==uuu) {
    printf("cand %u, deep %u, %u %u %u\n", u, deep, sg.deg[deep][u]-sg.pIdx[u], deep, edClique);
}
#endif
        if(sg.deg[deep][u]-sg.pIdx[u]+deep+edClique+1 < k) continue;
#ifdef DDEBUG
if(uu==uuu) {
    printf("cand %u, deep %u\n", u, deep);
}
#endif
        nxtC.clear();
        for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
            ui v = sg.pEdge[j];
#ifdef DDEBUG
if(uu==uuu) {
    printf("%u-%u:lv %u, deep%u\n", u, v, level[v], deep);
}
#endif
            // if(level[v] != deep) continue;
            nxtC.push_back(v);
            level[v] = deep + 1;
        }

        for(ui j = 0; j < nxtC.size(); j++) {
            ui v = nxtC[j];
            ui & ed = sg.deg[deep + 1][v];
            ed = sg.deg[deep][v];

            for(ui l = sg.pIdx[v]; l < ed; ) {
                ui w = sg.pEdge[l];
                if(level[w] == deep + 1) l++;
                else std::swap(sg.pEdge[l], sg.pEdge[--ed]);
            }
        }

        ui newEdClique = edClique;
        for(ui j = 0; j < newEdClique; ) {
            if(!g.connectHash(u, clique[j])) {
                std::swap(clique[j], clique[--newEdClique]);
            }
            else j++;
        }

        listClique(deep + 1, newEdClique);

        for(auto v : nxtC) level[v] = deep;
    }
}

kccPlusPivoter::kccPlusPivoter(Graph && g, ui k) :g(g), k(k) { 
    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answer = 0;

    computeC();

    sg.pIdx.resize(g.n);
    sg.pIdx2.resize(g.n);
    sg.pEdge.resize(g.m);

    sg.deg.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        sg.deg[i].resize(g.n);
    }

    nodes.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        nodes[i].resize(g.coreNumber);
    }

    level.resize(g.n);

    L.resize(g.n);
    for(ui i = 0; i < g.n; i++) {
        L[i].resize(g.coreNumber);
    }
    R.resize(g.n);
    for(ui i = 0; i < g.n; i++) {
        R[i].resize(g.coreNumber);
    }

    weight.resize(g.n);
    twoHopNodes.resize(g.coreNumber);
    inP.resize(g.n);

    degL.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        degL[i].resize(g.n);
    }
    degR.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        degR[i].resize(g.n);
    }

    printf("kccPlusPivoter2.cpp::kccPlusPivoter\n");
}
