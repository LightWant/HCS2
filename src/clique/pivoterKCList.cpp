#include "pivoter.h"

// #define DDEBUG
// #define BASELINE
#define BRANCHES

#ifdef BRANCHES
ui branches = 0;
#endif

#ifdef DDEBUG
ui uu = 0, uuu = 1;
#endif
std::vector<double> pivoter::run() {
    printf("pivoterKCList.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
#endif
        nodes[0].clear();
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            nodes[0].push_back(g.pEdge[i]);
            level[g.pEdge[i]] = 1;
        }
#ifdef DDEBUG
double p = answers[k];
#endif
        
        for(ui v : nodes[0]) {
            sg.pIdx[v] = sg.deg[0][v] = g.pIdx[v];

            for(ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(level[w] == 1) sg.pEdge[sg.deg[0][v]++] = w;
            }
        }
        for(ui v : nodes[0]) level[v] = 0;



        searchKCList(0, 0, 1);

#ifdef DDEBUG
printf("%u %.0f\n", u, answers[k] - p);
#endif
    }
#ifdef BASELINE
std::vector<ui> vc(g.n);
auto print = [&](uint32_t x) {
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
for(ui u = 0; u < g.n; u++) printf("%u %u\n", u, vc[u]);
#endif

#ifdef BRANCHES
printf("branches:%u\n", branches);
#endif

    return answers;
}

void pivoter::searchKCList(ui deep, ui p, ui h) {
#ifdef BRANCHES
branches++;
#endif

    if(h >= maxK) return;

    auto updateAns = [&]() {
        for(ui i = h; i < maxK && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
        }
    };

    std::vector<ui> & C = nodes[deep];
    ui edC = C.size();
#ifdef DDEBUG
if(uu == uuu) {
printf("    deep %u\nC:", deep);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("p %u, h %u\n", p, h);
}
#endif

    if(edC == 0) {
        updateAns(); return;
    }
    if(edC == 1) {
        p += 1; updateAns(); return;
    }
    if(edC == 2) {
        if(g.connectHash(C[0], C[1])) {
            p += 2; updateAns();
        }
        else {
            p += 1; updateAns();
            p -= 1; h += 1; updateAns();
        }
        return;
    }

    ui pivot = C[0], pivotDeg = 0, num = 0;
    for(ui i = 0; i < edC; i++) {
        ui v = C[i];
        ui tmp = 0;

        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(v, C[j])) tmp++;
        }

        if(tmp > pivotDeg) {
            pivot = v; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) {
            num++;
        }
    }
#ifdef DDEBUG
if(uu == uuu) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    if(pivotDeg+1 == edC && num == edC) {
        p += edC; updateAns(); return;
    }

    std::vector<ui> & nxtC = nodes[deep + 1];
    // C.changeTo(pivot, --edC);
    // ui newEdC = 0;
    nxtC.clear();
    for(ui i = 0; i < edC; i++) {
        if(C[i] < pivot && g.connectHash(pivot, C[i])) {
            nxtC.push_back(C[i]);
            level[C[i]] = deep + 1;
        }
    }
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

    searchKCList(deep + 1, p + 1, h);

    for(auto v : nxtC) level[v] = deep;
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
        if(g.connectHash(pivot, u)) {
#ifdef DDEBUG
if(uu == uuu) {
printf("    C1deep %u, cand %u, pivot %u\nC:", deep, u, pivot);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
            for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
                ui v = sg.pEdge[i];
                if(v == pivot || g.connectHash(pivot, v)) continue;
#ifdef DDEBUG
if(uu == uuu) {
printf("    C1deep %u, cand %u, godown with %u\nC:", deep, u, v);
}
#endif
                for(ui j = sg.pIdx[v]; j < sg.deg[deep][v]; j++) {
                    ui w = sg.pEdge[j]; 
                    ok[w] = true;
                }
                nxtC.clear();
                for(ui j = sg.pIdx[u]; j < sg.deg[deep][u]; j++) {
                    ui w = sg.pEdge[j];
                    if((w < v && g.connectHash(v, w) && g.connectHash(pivot, w)) || ok[w]) {
                        nxtC.push_back(w);
                        level[w] = deep + 1;
                    }
                }
                for(ui j = sg.pIdx[v]; j < sg.deg[deep][v]; j++) {
                    ui w = sg.pEdge[j]; ok[w] = false;
                }
                updateSG();
                searchKCList(deep + 1, p, h + 2);
                for(auto v:nxtC) level[v] = deep;
            }
        }
        else {
#ifdef DDEBUG
if(uu == uuu) {
printf("    C2deep %u, cand %u, pivot %u\nC:", deep, u, pivot);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
            nxtC.clear();
            for(ui i = sg.pIdx[u]; i < sg.deg[deep][u]; i++) {
                ui v = sg.pEdge[i];
                // if(v == pivot) continue;
                nxtC.push_back(v);
                level[v] = deep + 1;
            }
            updateSG();
            searchKCList(deep + 1, p, h + 1);
            for(auto v:nxtC) level[v] = deep;
        }
    }
}

pivoter::pivoter(Graph && g, ui k) :g(g), k(k) { 
    C.resize(g.n);

    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answers.resize(maxK, 0.0);

    computeC();

    initBuffer(g.coreNumber * maxK);

    sg.pIdx.resize(g.n);
    sg.pIdx2.resize(g.n);
    sg.pEdge.resize(g.m);
    sg.deg.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        sg.deg[i].resize(g.n);
    }

    nodes.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        nodes[i].resize(g.coreNumber);
    }

    level.resize(g.n);
    ok.resize(g.n, false);

    printf("pivoter.h\n");
}

