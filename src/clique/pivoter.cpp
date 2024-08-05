#include "pivoter.h"
#include <chrono>
#include <ctime>
#define DDEBUG
// #define BASELINE
#define BRANCHES

#ifdef BRANCHES
ui branches = 0;
#endif

#ifdef DDEBUG
ui uu = 0, d = 0, uuu=189;
#endif

std::vector<double> pivoter::run() {
    printf("pivoter.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
#endif
// auto s1 = std::chrono::steady_clock::now();
// double s1 = clock();
if(u==189) {
// // if(u==1696186) {//ski
// for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
//     ui v = g.pEdge[i];
//     for(ui j = i + 1; j < g.pIdx[u + 1]; j++) {
//         ui w = g.pEdge[j];
//         if(g.connectHash(v, w)) printf("%u %u\n", v, w);
//     }
// }
}
// else continue;
        ui edC = 0;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], edC++);
        }
#ifdef DDEBUG
double p = answers[k];
printf("upt %u %.0f\n", u, answers[k]);
#endif
// double p = answers[k];
double p1 = answers[k-1];
        search(edC, 0, 1);
#ifdef DDEBUG
printf("ut %u %.0f ak%.0f\n", u, answers[k] - p, answers[k]);
#endif
if(u==189) {
    for(ui i = 0; i < edC; i++) {
        printf("%u ", C[i]);
    }printf("\n");
    printf("%u edC, %u maxK\n", edC, maxK);
    bool isClique = true;
    for(ui i = 0; i < edC; i++) {
        for(ui j = i + 1; j < edC; j++) {
            if(!g.connectHash(C[i], C[j])) isClique = false;
        }
    }
    if(isClique) printf("isCLique\n");
}
// auto s2 = std::chrono::steady_clock::now();
// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
// if(duration.count() >= 10) {
//     printf("%u %.0f\n", u, (double)duration.count());
// }
// double s2 = clock();
// if((s2 - s1) / CLOCKS_PER_SEC > 0.01) {
//     printf("%u\n", u);fflush(stdout);
// }
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

void pivoter::search(ui edC, ui p, ui h) {
#ifdef BRANCHES
branches++;
#endif
#ifdef DDEBUG
d++;
if(uu == uuu) {
printf("    deep %u\nC:", d);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("p %u, h %u\n", p, h);
}
#endif
    if(h >= maxK) return;
    auto updateAns = [&]() {
#ifdef DDEBUG
if(uu == uuu) {
printf("updataAns p%u h%u maxK%u, C[p][29] %.0f\n", p, h, maxK, CN[p][29]);
}
#endif
        for(ui i = h; i < maxK && i <= p + h; i++) {
#ifdef DDEBUG
if(uu == uuu) {
printf("ans[%u] %.0f+=CN[%u][%u] %.0f\n", i, answers[i], p, i-h, CN[p][i-h]);
}
#endif
            answers[i] += CN[p][i-h];
#ifdef DDEBUG
if(uu == uuu) {
printf("ans[%u] %.0f\n", i, answers[i]);
}
#endif
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

    ui pivot = C[0], pivotDeg = 0, num = 0;
    for(ui i = 0; i < edC; i++) {
        ui tmp = 0;

        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(C[i], C[j])) tmp++;
        }

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
#ifdef DDEBUG
d--;
if(uu == uuu) {
printf("isCLique add %u edC to p\n", edC);
}
#endif
        p += edC; updateAns(); return;
    }
// if(uu >= de) {
//     printf("p %u pd %u\n", pivot, pivotDeg); fflush(stdout);
// }
    C.changeTo(pivot, --edC);
    ui newEdC = 0;
    for(ui i = 0; i < edC; i++) {
        if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++); 
    }

    search(newEdC, p + 1, h);
#ifdef DDEBUG
if(uu == uuu) {
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
if(uu == uuu) {
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
    C.resize(g.n);

    // sg.pIdx.resize(g.n);
    // sg.pIdx2.resize(g.n);
    // sg.pEdge.resize(g.m);

    // sg.deg.resize(g.n);
    // sg.deg[0].resize(g.n);

    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answers.resize(maxK, 0.0);

    computeC();

    initBuffer(g.coreNumber * maxK);

    printf("pivoter.h\n");
}
