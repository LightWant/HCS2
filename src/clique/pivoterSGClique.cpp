#include "pivoter.h"
#include <cassert>

// #define DDEBUG
// #define DEBUG
// #define BASELINE
// #define BRANCHES

#ifdef BRANCHES
ui branches = 0;
ui PNODES = 0, HNODES = 0;
ui cntPQ[500][500];
#endif

#ifdef DDEBUG
ui uu = 0, uuu = 118;
#endif
std::vector<double> pivoter::run() {
    printf("pivoterSGClique.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
// g.print();
#endif

    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;

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
double p = answers[k];
#endif
// double p = answers[k];
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
        for(ui i = 0; i < clique.size(); i++) {
            for(ui j = i + 1; j < clique.size(); j++) {
                assert(g.connectHash(clique[i], clique[j]));
            }
        }
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
        // for(ui i = 2; i+1 < maxK && i < clique.size(); i++) {
        //     answers[i+1] += CN[clique.size()][i];
        // }
#ifdef DDEBUG
sc4.clear();
#endif
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
        Ps.clear(); Hs.clear();
        Hs.push_back(u);
        searchSGClique(deep, 0, 1, clique.size());

        for(ui v : nodes[deep]) level[v] = 0;
// printf("ut %u %.0f\n", u, answers[k] - p);
#ifdef DDEBUG
printf("%u %.0f\n", u, answers[k] - p);
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
printf("hNodes:%u\npNodes:%u\n", HNODES, PNODES);
for(ui i = 0; i < g.coreNumber; i++) {
    for(ui j = 0; j < g.coreNumber; j++) {
        printf("%u ", cntPQ[i][j]);
    }
    printf("\n");
}
#endif
    return answers;
}

void pivoter::searchSGClique(ui deep, ui p, ui h, ui edClique) {
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
}
#endif
    if(h >= maxK) return;

    auto updateAns = [&]() {
#ifdef BRANCHES
cntPQ[p][edClique]++;
return;
#endif

#ifdef DDEBUG
if(uu == uuu) {
printf("\n In updateAns\n");
printf("    deep %u\nC:", deep);
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
}
#endif

//assert clique
// fflush(stdout);
// for(ui i = 1; i < Hs.size(); i++) {
//     for(ui j = 0; j < edClique; j++) {
//         assert(g.connectHash(Hs[i], clique[j]));
//     }
//     bool f = true;
//     for(ui j = edClique; j < clique.size(); j++) {
//         if(!g.connectHash(Hs[i], clique[j])) {
//             f = false; break;
//         }
//     }
//     if(clique.size() == edClique) f = false;
//     assert(!f);
// }  
#ifdef DDEBUG  
auto buildSC4 = [&](ui p, std::vector<ui> & Ps) {
for(ui m = (1<<p)-1; m > 0; m--) {
    ui sz = 0;
    for(ui i = 0; i < p; i++) {
        if((1<<i) & m) sz++;
    }
    if(sz+h != k) continue;
    c4 tmp;
    for(ui i = 0; i < h; i++) tmp.v[i] = Hs[i];
    ui j = h;
    for(ui i = 0; i < p; i++) {
        if((1<<i) & m) {
            tmp.v[j++] = Ps[i];
        }
    }
    assert(j == k);
    tmp.sort();
    if(!tmp.isClique(&g)) {
        printf("noclique_type2:");
        for(ui i = 0; i < k; i++) printf("%u ", tmp.v[i]); printf("\n");
    }
    if(sc4.find(tmp) != sc4.end()) {
        printf("multi1 clique:");
        for(ui i = 0; i < k; i++) printf("%u ", tmp.v[i]); printf("\n");
        fflush(stdout);
    }
    else sc4.insert(tmp);
}
};
#endif
        for(ui i = h; i < maxK && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
#ifdef DDEBUG
if(uu == uuu) {
if(i==k) {
printf("addp bi to ans:%u %.0f\n",
     i, CN[p][i-h]);

if(k==4)buildSC4(p, Ps);
}
}
#endif
        }
        for(ui i = 1; i+h < maxK && i <= edClique; i++) {
            answers[i + h] += CN[edClique][i];
#ifdef DDEBUG
if(uu == uuu) {
if(i+h==k) {
printf("addp clique to ans:i%u, edCique %u %.0f\n",
        i, edClique, CN[edClique][i]);

if(k==4)
buildSC4(edClique, clique);
}
}
#endif
        }

        if(p == 0 || edClique == 0) return;
        //build the bipartite graph
        ui pL = 0, pR = 0;
        for(ui i = 0; i < p; i++) {
            bool connected = false;
            for(ui j = 0; j < edClique; j++) {
                if(g.connectHash(Ps[i], clique[j])) {
                    connected = true;
                    break;
                }
            }
            if(connected) candL.changeTo(Ps[i], pL++);
        }
        for(ui i = 0; i < edClique; i++) {
            bool connected = false;
            for(ui j = 0; j < p; j++) {
                if(g.connectHash(clique[i], Ps[j])) {
                    connected = true;
// #ifdef DDEBUG
// if(uu == uuu) {
// printf("connected %u %u\n", clique[i], Ps[j]);
// }
// #endif
                    break;
                }
            }
            if(connected) candR.changeTo(clique[i], pR++);
        }

        //choose a=k-h nodes in total
        //(1,a-1)-biclique, (2,a-2)-biclique
        //(3,a-3)-biclique.....(a-1,1)-biclique
        for(ui i = 2; i < maxK && i <= pL + pR; i++) 
            ansBiclique[i] = 0;
        limitOfPQ = maxK-1;
        // pivotCount(0, pL, pR, {0,0,0,0});//biclique counting

        for(ui i = 2; i + h < maxK && i <= pL + pR; i++) {
#ifdef DDEBUG
double pp = answers[i+h];
#endif
            answers[i + h] += ansBiclique[i];
#ifdef DDEBUG
if(uu == uuu) {
ui sumDeg = 0, k3 = 0, k4 = 0;
assert(answers[i+h] - ansBiclique[i] == pp);

if(i+h==k) {
    printf("add2 ansbi to ans:i+h%u i%u %u\n", i + h, i, ansBiclique[i]);
    // double cnti = getBicliqueSizeI(i);
}

}
#endif


        }
#ifdef DDEBUG
if(uu == uuu) {
printf("            ans[k] %.0f\n", answers[k]);
}
#endif

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
            p += 1; updateAns();
            Ps.pop_back();
            Hs.push_back(C[1]);
            p -= 1; h += 1; 
            for(ui j = 0; j < edClique; ) {
                if(!g.connectHash(C[1], clique[j])) {
                    std::swap(clique[j], clique[--edClique]);
                }
                else j++;
            }
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
        for(ui i = 0; i < edC; i++) Ps.push_back(C[i]);
#ifdef DDEBUG
if(uu == uuu) {
printf("cand is a clique\n");
}
#endif
        updateAns(); 
        for(ui i = 0; i < edC; i++) Ps.pop_back();
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

#ifdef BRANCHES
PNODES++;
#endif
    Ps.push_back(pivot);
    searchSGClique(deep + 1, p + 1, h, edClique);
    Ps.pop_back();
    for(auto v:nxtC) level[v] = deep;

#ifdef DDEBUG
if(uu == uuu) {
printf("    redeep %u, \nC:", deep);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
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
            if(level[w] == deep - 1) continue;
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
#ifdef BRANCHES
HNODES++;
#endif
        searchSGClique(deep + 1, p, h + 1, newEdClique);
        Hs.pop_back();

        for(auto v : nxtC) level[v] = deep;
    }
}

pivoter::pivoter(Graph && g, ui k) :g(g), k(k) { 
    maxSize = g.coreNumber + 5;
    maxK = k + 1;
    answers.resize(maxK, 0.0);
    ansBiclique.resize(maxK, 0.0);

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

    tmpNodesL.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        tmpNodesL[i].resize(g.coreNumber);
    }
    tmpNodesR.resize(maxSize);
    for(ui i = 0; i < maxSize; i++) {
        tmpNodesR[i].resize(g.coreNumber);
    }

    initBuffer(g.coreNumber * maxK * maxK);

    level.resize(g.n);
   
    candL.resize(g.n);
    candR.resize(g.n);

    printf("pivoter.h\n");
}

void pivoter::pivotCount(int l, int pL, int pR, treePath t) {
#ifdef BRANCHES
branches++;
#endif
#ifdef DEBUG
if(uu == uuu) {
printf("\nl:%d, pL:%d, pR:%d, %d :%d :%d :%d\n",
    l, pL, pR, t.p1, t.h1, t.p2, t.h2);
printf("candL:");
for(int i = 0; i < pL; i++) {
printf("%u ", candL[i]);
} printf("\n");

printf("candR:");
for(int i = 0; i < pR; i++) {
    printf("%u ", candR[i]);
}
printf("\n\n");
}
#endif
#ifdef DEBUG
int x = 3, y = 1;
#endif
// if(l > 10) return;

    auto updateAnsBiclique = [&]() {
#ifdef DEBUG
if(uu == uuu) {
printf("%d :%d :%d :%d\n", t.p1, t.h1, t.p2, t.h2);
}
#endif
        if(t.p1 + t.h1 == 0 || t.p2 + t.h2 == 0) return;
        for(int ll = 0; ll <= t.p1; ll++) {
            for(int r = 0; r <= t.p2; r++) {
                if(ll+t.h1 + r+t.h2 > limitOfPQ) break;
                if(ll + t.h1 == 0 || r + t.h2 == 0) continue;
                // if(!(ll + t.h1 == 2 && r + t.h2 == 2)) continue;
#ifdef DEBUG
if(uu == uuu) {
if(ll + t.h1 + r + t.h2 == limitOfPQ) {
    printf("    add asbiclique k-1, %.0f\n", CN[t.p1][ll] * CN[t.p2][r]);
}
}
#endif
                ansBiclique[ll + t.h1 + r + t.h2] 
                    += CN[t.p1][ll] * CN[t.p2][r];
            }
        }
    };
    if(t.h1 + t.h2 >= maxK) return;
    if(pL == 0 && pR == 0) {
        // t.p1 += pL;
        // t.p2 += pR;
#ifdef DEBUG
if(uu == uuu) {
printf("adding p1 %d h1 %d p2 %d h2 %d\n", t.p1, t.h1, t.p2, t.h2);
}
#endif
        updateAnsBiclique();
    
        return ;
    }

    int pivotL = 0;
    int maxDeg = 0;
    for(int i = 0; i < pL; i++) {
        int deg = 0;

        uint32_t u = candL[i];

        for(int j = 0; j < pR; j++) {
            if(g.connectHash(u, candR[j])) {
                deg++;
            }
        }

        if(deg > maxDeg) {
            maxDeg = deg;
            pivotL = u;
        }
    }
#ifdef DEBUG
if(uu == uuu) {
printf("pivotL:%d, %d\n", pivotL, maxDeg);
}
#endif

    if(maxDeg == 0) {
        pivotCount(l + 1, 0, 0, {t.p1 + pL, t.h1, t.p2, t.h2});

        
        for(int i = 1; i <= pR; i++) {
            // pivotCount(l + 1, 0, 0, {t.p1, t.h1, t.p2 + pR, t.h2 + i})*C[pR][i];
            t.h2 += i;
#ifdef DEBUG
if(uu == uuu) {
printf("xadding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
}
#endif
            // updateAnsBiclique();
#ifdef DEBUG
if(uu == uuu) {
printf("%d :%d :%d :%d\n", t.p1, t.h1, t.p2, t.h2);
}
#endif
            if(t.p1 + t.h1 == 0 || t.p2 + t.h2 == 0) return;
            for(int ll = 0; ll <= t.p1; ll++) {
                for(int r = 0; r <= t.p2; r++) {
                    if(ll+t.h1 + r+t.h2 >= limitOfPQ) break;
                    if(ll + t.h1 == 0 || r + t.h2 == 0) continue;
                    // if(!(ll + t.h1 == 2 && r + t.h2 == 2)) continue;
#ifdef DEBUG
if(uu == uuu) {
if(ll + t.h1 + r + t.h2 == limitOfPQ) {
    printf("    add asbiclique k-1*C[pR][i], %.0f\n", 
        CN[t.p1][ll] * CN[t.p2][r]* CN[pR][i]);
}
}
#endif
                    ansBiclique[ll + t.h1 + r + t.h2] 
                        += CN[t.p1][ll] * CN[t.p2][r] * CN[pR][i];
                }
            }

            t.h2 -= i;
        }
    
        return;
    }

    int pRR = 0;
    for(int j = 0; j < pR; j++) {
        if(g.connectHash(pivotL, candR[j])) {
            candR.changeToByPos(j, pRR++);
        }
    }
    int pivotR = 0;
    maxDeg = 0;
    for(int i = 0; i < pRR; i++) {
        int deg = 0;
        int v = candR[i];

        for(int j = 0; j < pL; j++) {
            if(g.connectHash(v, candL[j])) {
                deg++;
            }
        }

        if(deg > maxDeg) {
            maxDeg = deg;
            pivotR = v;
        }
    }
    if(pRR > 0) candR.changeTo(pivotR, --pRR);
    int pLL = 0;
    for(int j = 0; j < pL; j++) {
        if(g.connectHash(pivotR, candL[j])) {
            candL.changeToByPos(j, pLL++);
        }
    }
    if(pLL > 0) candL.changeTo(pivotL, --pLL);

#ifdef DEBUG
if(uu == uuu) {
printf("choose p %u-%u\n", pivotL, pivotR);
}
#endif

    pivotCount(l + 1, pLL, pRR, 
        treePath{t.p1 + 1, t.h1, t.p2 + 1, t.h2});

    // if((int)tmpNodesL[l].capacity() < pL - pLL) tmpNodesL[l].resize((pL - pLL)*1.5);
    candL.copy(tmpNodesL[l].data(), pL, pLL + 1);
    // if((int)tmpNodesR[l].capacity() < pR - pRR) tmpNodesR[l].resize((pR - pRR)*1.5);
    candR.copy(tmpNodesR[l].data(), pR, pRR + 1);

    int tmpLSize = pL - pLL - 1;
    int tmpRSize = pR - pRR - 1;
#ifdef DEBUG
if(uu == uuu) {
printf("l%u, tmpLRSize:%d,%d, pLR:%d,%d\n", l, tmpLSize, tmpRSize, pL, pR);
}
#endif

    for(int i = 0; i < tmpLSize; i++) {
        uint32_t u = tmpNodesL[l][i];
        candL.changeTo(u, --pL);
        pRR = 0;

        // for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
        //     uint32_t v = g->e1[j];
        //     if(candR.idx(v) < (uint32_t)pR) {
        //         candR.changeTo(v, pRR++);
        //     }
        // }
        for(int j = 0; j < pR; j++) {
            if(g.connectHash(u, candR[j])) {
                candR.changeToByPos(j, pRR++);
            }
        }
#ifdef DEBUG
if(uu == uuu) {
printf("testing h %u, pRR %d\n", u, pRR);
}
#endif
        t.p1 += pL;
        t.h1 += 1;
        // for(int j = 1; j <= tmpRSize - i; j++) {
            // t.h2 += j;
#ifdef DEBUG
if(uu == uuu) {
printf("x1adding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
}
#endif

        updateAnsBiclique();

            // t.h2 -= j;
        // }
        t.p1 -= pL;
        t.h1 -= 1;

        ui totalSz = pR;
        ui * cand = allocMem(totalSz);
        ui candSize = 0;
        for(int i = 0; i < pR; i++) {
            if(g.connectHash(u, candR[i])) {
                cand[candSize++] = candR[i];
            }
        }
        for(ui i = 0; i < candSize; i++) {
            ui v = cand[i];
            pLL = 0;
            for(int k = 0; k < pL; k++) {
                if(g.connectHash(v, candL[k])) {
                    candL.changeToByPos(k, pLL++);
                }
            }
            candR.changeTo(v, --pRR);
#ifdef DEBUG
if(uu == uuu) {
printf("choose h %u-%u\n", u, v);
}
#endif

            pivotCount(l + 1, pLL, pRR, 
                treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
        }
        freeMem(totalSz);
    }


    for(int i = 0; i < tmpRSize; i++) {
        uint32_t v = tmpNodesR[l][i];
        candR.changeTo(v, --pR);
        pLL = 0;

        // for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
        //     uint32_t u = g->e2[j];
        //     if(candL.idx(u) < (uint32_t)pL) {
        //         candL.changeTo(u, pLL++);
        //     }
        // }
        for(int j = 0; j < pL; j++) {
            if(g.connectHash(v, candL[j])) {
                candL.changeToByPos(j, pLL++);
            }
        }

        t.p2 += pR;//pivotR + [i+1,tmpRSize-1]
        t.h2 += 1;
        // for(int j = 1; j <= tmpRSize - i; j++) {
            // t.h2 += j;
#ifdef DEBUG
if(uu == uuu) {
printf("x2adding %d %d %d %d, %d/%d\n", t.p1, t.h1, t.p2, t.h2, i, pR);
}
#endif

        // for(int l = 0; l <= t.p1 && l < maxD2; l++) {
        //     for(int r = 0; r <= t.p2 && r < maxD2; r++) {
        //         if(r + t.h2 >= (int)ansAll[l + t.h1].size()) {
        //             ansAll[l + t.h1].resize(r + t.h2 + 5);
        //         }
        //         ansAll[l + t.h1][r + t.h2] 
        //             += C[t.p1][l] * C[t.p2][r];
        //     }
        // }
        updateAnsBiclique();
            // t.h2 -= j;
        // }
        t.p2 -= pR;
        t.h2 -= 1;

        ui totalSz = pL;
        ui * cand = allocMem(totalSz);
        ui candSize = 0;
        for(int i = 0; i < pL; i++) {
            if(g.connectHash(v, candL[i]) && g.connectHash(pivotR, candL[i])) {
                cand[candSize++] = candL[i];
            }
        }
        for(ui i = 0; i < candSize; i++) {
            ui u = cand[i];
            pRR = 0;
            for(int k = 0; k < pR; k++) {
                if(g.connectHash(u, candR[k])) {
                    candR.changeToByPos(k, pRR++);
                }
            }
            candL.changeTo(u, --pLL);
#ifdef DEBUG
if(uu == uuu) {
printf("choose h %u-%u\n", u, v);
}
#endif
            pivotCount(l + 1, pLL, pRR, 
                treePath{t.p1, t.h1 + 1, t.p2, t.h2 + 1});
        }
        freeMem(totalSz);
    }
}


