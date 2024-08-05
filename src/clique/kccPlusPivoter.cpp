#include "kccPlusPivoter.h"
#include <cassert>

// #define DDEBUG
// #define DEBUG
// #define BASELINE
#define BRANCHES

#ifdef BRANCHES
ui branches = 0;
ui cntPQ[500][500];
#endif

#ifdef DDEBUG
ui uu = 0, uuu = 0;
#endif
double kccPlusPivoter::run() {
    printf("kccPlusPivoter.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
if(uu=uuu)
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
        // for(ui i = 2; i+1 < maxK && i < clique.size(); i++) {
        //     answers[i+1] += CN[clique.size()][i];
        // }

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

        // for(ui v : nodes[deep]) level[v] = deep;
        // listClique(deep, clique.size());

        for(ui v : nodes[deep]) level[v] = 0;
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
for(ui i = 0; i < g.coreNumber; i++) {
    for(ui j = 0; j < g.coreNumber; j++) {
        printf("%u ", cntPQ[i][j]);
    }
    printf("\n");
}
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
    if(p+h+nodes[deep].size()+2 < k) return;
    if(nodes[deep].size() + p + h + edClique < k) return;
    
    //choose 0,1,2 vertices from clique
    auto updateAns = [&]() {
#ifdef BRANCHES
cntPQ[p][edClique]++;
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
        for(ui i = 0; i < edClique; i++) L[i].clear();
        for(ui i = 0; i < p; i++) R[i].clear();
        for(ui i = 0; i < edClique; i++) {
#ifdef DDEBUG
if(uu == uuu) {
    printf("%u:", clique[i]);
}
#endif
            for(ui j = 0; j < p; j++) {
                if(g.connectHash(clique[i], Ps[j])) {
                    L[i].push_back(j);
                    R[j].push_back(i);
#ifdef DDEBUG
if(uu == uuu) {
    printf("%u ", Ps[j]);
}
#endif
                }
            }
        }
#ifdef DDEBUG
if(uu == uuu) {
    printf("\n");
}
#endif
        //choose 0 vertex from clique
        if(p >= k - h) answer += CN[p][k - h];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans0 += C(%u,%u)\n", p, k-h);
}
#endif
        //choose 1 vertex from clique
        for(ui i = 0; i < edClique; i++) {
            if(L[i].size() >= k-h-1) {
                answer += CN[L[i].size()][k-h-1];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans1 += C(%u,%u), clique node %u\n", L[i].size(), k-h-1, clique[i]);
}
#endif
            }
        }

        // return;

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
        for(ui i = 0; i < edClique; i++) weight[i] = 0;
        twoHopNodes.clear();
        for(ui i = 0; i < edClique; i++) {
            for(auto j : L[i]) {
                if(R[j].size())
                for(int l = R[j].size() - 1; l >= 0; l--) {
                    if(R[j][l] <= i) break;
                    weight[R[j][l]]++;
                    if(!vis[R[j][l]]) {
                        vis[R[j][l]] = true;
                        twoHopNodes.push_back(R[j][l]);
                    }
                }
            }

            for(auto v : twoHopNodes) {
                vis[v] = false;
                if(weight[v] >= k-h-2) {
                    answer += CN[weight[v]][k-h-2];
#ifdef DDEBUG
if(uu == uuu) {
    printf("ans2 += C(%u,%u), clique node %u-%u\n", weight[v], k-h-2, clique[i], clique[v]);
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
        searchSGClique(deep + 1, p, h + 1, newEdClique);
        Hs.pop_back();

        for(auto v : nxtC) level[v] = deep;
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
    if(deep == k-2) {
        if(edClique >= 2)
            answer += CN[edClique][2];
        return;
    }

    std::vector<ui> & C = nodes[deep];

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
        if(sg.deg[deep][u]-sg.pIdx[u]+deep+edClique+1 < k) continue;
        level[u] = deep-1;
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
            if(level[v] != deep) continue;
            nxtC.push_back(v);
            level[v] = deep + 1;
        }

        for(ui j = 0; j < nxtC.size(); j++) {
            ui v = nxtC[j];
            ui & ed = sg.deg[deep + 1][v];
            ed = sg.deg[deep][v];

            for(ui l = sg.deg[deep][v]; l < ed; ) {
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
    maxSize = g.coreNumber + 5;
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

    L.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        L.resize(g.coreNumber);
    }
    R.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        R.resize(g.coreNumber);
    }

    weight.resize(g.coreNumber);
    twoHopNodes.resize(g.coreNumber);
    vis.resize(g.coreNumber);

    printf("kccPlusPivoter.cpp::kccPlusPivoter\n");
}
