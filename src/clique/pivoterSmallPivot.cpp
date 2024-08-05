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
ui uu = 0, d = 0, uuu = 2;
#endif
std::vector<double> pivoter::run() {
    printf("pivoterSmallPivot.cpp::run\n");
    g.initHash();
#ifdef DDEBUG
g.print();
#endif
    for(ui u = 0; u < g.n; u++) {
#ifdef DDEBUG
uu = u;
#endif

        ui edC = 0;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], edC++);
        }
#ifdef DDEBUG
double p = answers[k];
#endif
        search(0, edC, 0, 1);
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
printf("branches:%u\n", branches); fflush(stdout);
#endif
    return answers;
}

void pivoter::search(ui stC, ui edC, ui p, ui h) {
#ifdef BRANCHES
branches++;
#endif
#ifdef DDEBUG
d++;
if(uu == uuu) {
printf("    deep %u\nC:", d);
for(ui i = stC; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("p %u, h %u\n", p, h);
}

#ifdef DDEBUG
if(edC-stC <= 2 || h >= maxK) d--;
#endif
#endif
    if(h >= maxK) return;
    auto updateAns = [&]() {
        for(ui i = h; i < maxK && i <= p + h; i++) {
            answers[i] += CN[p][i-h];
        }
    };

    if(edC == stC) {
        updateAns(); return;
    }
    if(edC == 1 + stC) {
        p += 1; updateAns(); return;
    }
    if(edC == 2 + stC) {
        if(g.connectHash(C[stC], C[stC + 1])) {
            p += 2; updateAns();
        }
        else {
            p += 1; updateAns();
            p -= 1; h += 1; updateAns();
        }
        return;
    }

    while(edC > stC) {
        ui pivot = C[stC], pivotDeg = edC - stC - 1, num = 1;
        for(ui i = stC; i < edC; i++) {
            ui tmp = 0;

            for(ui j = stC; j < edC; j++) {
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

        if(pivotDeg == edC - stC - 1 && num == edC - stC) {
#ifdef DDEBUG
d--;
#endif
            p += edC - stC; updateAns(); return;
        }

        C.changeTo(pivot, --edC);
        ui newEdC = stC;
        for(ui i = stC; i < edC; i++) {
            if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++); 
        }

        search(stC, newEdC, p + 1, h);
#ifdef DDEBUG
if(uu == uuu) {
printf("    redeep %u, \nC:", d);
for(ui i = 0; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
        ui candSize = newEdC - stC;
        ui * cand = allocMem(candSize);
        memcpy(cand, C.begin() + stC, sizeof(ui) * candSize);
    
        for(ui i = 0; i < candSize; i++) {
            ui v = cand[i];
            C.changeTo(v, --edC);

            newEdC = stC;
            for(ui j = stC; j < edC; j++) {
                if(g.connectHash(v, C[j])) C.changeToByPos(j, newEdC++); 
            }
#ifdef DDEBUG
if(uu == uuu) {
printf("    deep %u, i %u, candv %u\nC:", d, i, v);
for(ui i = stC; i < edC; i++) {
printf("%u ", C[i]);
}printf("\n");
printf("newC:");
for(ui i = stC; i < newEdC; i++) {
printf("%u ", C[i]);
}printf("\n");
}
#endif
            ui candSize2 = newEdC - stC;
            ui * cand2 = allocMem(candSize2);
            ui nonNeiOfPivot = 0;
            for(ui j = stC; j < newEdC; j++) {
                ui w = C[j];
                if(g.connectHash(pivot, w)) continue;
                cand2[nonNeiOfPivot++] = w;
            }
            ui newStC = stC;
            for(ui j = 0; j < nonNeiOfPivot; j++) {
                ui w = cand2[j];
                C.changeTo(w, newStC++);
#ifdef DDEBUG
if(uu == uuu) {
printf("    GoDown deep %u, candv %u-%u\nC:", d, v, w);
}
#endif
                ui newEdC2 = newStC;
                for(ui l = newStC; l < newEdC; l++) {
                    if(g.connectHash(w, C[l])) 
                        C.changeToByPos(l, newEdC2++);
                }

                search(newStC, newEdC2, p, h + 2);
            }
            freeMem(candSize2);
        }

        freeMem(candSize);
    }
#ifdef DDEBUG
d--;
#endif
}

pivoter::pivoter(Graph && g, ui k) :g(g), k(k) { 
    C.resize(g.n);

    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answers.resize(maxK, 0.0);

    computeC();

    initBuffer(g.coreNumber * maxK);

    printf("pivoter.h\n");
}
