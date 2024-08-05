#include "pivoter.h"
#include <chrono>
#include <ctime>
// #define DDEBUG
// #define BASELINE
#define BRANCHES

#ifdef BRANCHES
ui branches = 0;
#endif

#ifdef DDEBUG
ui uu = 0, d = 0;
#endif
std::vector<double> pivoter::run() {
    printf("pivoterClique.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
#endif
        constexpr ui deep = 1;
        nodes[deep].clear();
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            nodes[deep].push_back(g.pEdge[i]);
            level[g.pEdge[i]] = deep;
        }
#ifdef DDEBUG
double p = answers[k];
#endif
        if(nodes[deep].size() == 0) continue;

        ui maxDeg = 0, maxV = C[0];
        findTheMax(0, edC, maxDeg, maxV);

        searchSGClique(deep, 0, stC, stC, edC, 0, 1);
#ifdef DDEBUG
printf("%u %.0f\n", u, answers[k] - p);
#endif
    }
#ifdef BASELINE
auto print = [&](uint32_t x) {
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");
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

#endif

#ifdef BRANCHES
printf("branches:%u\n", branches); fflush(stdout);
#endif
    return answers;
}

void pivoter::findTheMax(ui st, ui ed, ui & maxDeg, ui & maxV) {
    for(ui i = st; i < ed; i++) {
        ui deg = 0;
        for(ui j = st; j < ed; j++) {
            if(g.connectHash(C[i], C[j])) deg++;
        }
        if(deg > maxDeg) {
            maxDeg = deg;
            maxV = C[i];
            if(deg + 1 == ed - st) break;
        }
    } 
}
ui pivoter::expandClique(ui stC, ui edC) {
//stC-edC is a clique, expand it
    ui newStC = stC;
    for(ui i = 0; i < newStC; ) {
        bool isCommonNei = true;
        for(ui j = stC; j < edC; j++) {
            if(!g.connectHash(C[i], C[j])) {
                isCommonNei =  false;
                break;
            }
        }
        if(isCommonNei) C.changeToByPos(i, --newStC);
        else i++;
    }

    while(newStC < stC) {
        ui maxDeg = 0, maxV = C[newStC];
        findTheMax(newStC, stC, maxDeg, maxV);
        C.changeTo(maxV, --stC);
        ui tmp = stC;
        for(ui i = newStC; i < tmp; ) {
            if(g.connectHash(maxV, C[i])) C.changeToByPos(i, --tmp);
            else i++;
        }
        newStC = tmp;
    }

    return stC;
}

void pivoter::searchClique(ui st, ui ed, ui stC, ui edC, ui p, ui h) {
#ifdef BRANCHES
branches++;
#endif
#ifdef DDEBUG
d++;
if(uu == 2) {
printf("    deep %u\nC:", d);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("p %u, h %u\n", p, h);
}
#endif
    if(h >= maxK) return;
    auto updateAns = [&]() {
        for(ui i = h; i < maxK && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
        }
    };
#ifdef DDEBUG
if(edC <= 2) d--;
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

    ui nstC = expandClique(stC, edC);

    ui pivot = C[0], pivotDeg = 0;
    findTheMax(st, stC, pivotDeg, pivot);
    C.changeTo(pivot, st++);

    ui newSt = stC;
    if(stC > 0)
    for(ui i = stC - 1; i >= st; i--) {
        if(g.connectHash(pivot, C[i])) C.changeToByPos(i, --newSt); 
    }
    ui newEdC = stC;
    for(ui i = stC; i < edC; i++) {
        if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++); 
    }

#ifdef DDEBUG
if(uu == 2) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    searchClique(newSt, stC, newEdC, p + 1, h);
#ifdef DDEBUG
if(uu == 2) {
printf("    redeep %u, \nC:", d);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
    ui candSize = edC - pivotDeg;
    ui * cand = allocMem(candSize);
    memcpy(cand, C.begin() + pivotDeg, sizeof(ui) * candSize);

    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];
        C.changeTo(v, --edC);
#ifdef DDEBUG
if(uu == 2) {
printf("    deep %u, i %u, candv %u\nC:", d, i, v);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
        ui newEdC = 0;
        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(v, C[j])) C.changeToByPos(j, newEdC++); 
        }

        search(newEdC, p, h + 1);
    }

    freeMem(candSize);
#ifdef DDEBUG
d--;
#endif
}

pivoter::pivoter(Graph && g, ui k) :g(g), k(k) { 
    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answers.resize(maxK, 0.0);

    computeC();

    initBuffer(g.coreNumber * maxK);

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
    // ok.resize(g.n, false);

    printf("pivoter.h\n");
}
