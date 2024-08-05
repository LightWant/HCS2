#include "pivoter.h"

// #define DDEBUG
// #define BASELINE

#ifdef DDEBUG
ui uu = 0, d = 0;
#endif
std::vector<double> pivoter::run() {
    printf("pivoter2.cpp::run\n");
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
double p = answers[4];
#endif
        // if(edC < 20) {
        //     search(edC, 0, 1);
        //     continue;
        // }
        
        for(ui i = 0; i < edC; i++) {
            ui v = C[i];
            sg.pIdx[v] = sg.deg[0][v] = g.pIdx[v];
            for(ui i = 0; i < edC; i++) {
                if(g.connectHash(v, C[i])) sg.pEdge[sg.deg[0][v]++] = C[i];
            }
        }
        searchSG(0, edC, 0, 1);

#ifdef DDEBUG
printf("%u %.0f\n", u, answers[4] - p);
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
    return answers;
}

void pivoter::searchSG(ui deep, ui edC, ui p, ui h) {
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

    ui pivot = C[0], pivotDeg = 0, num = 0;
    for(ui i = 0; i < edC; i++) {
        // ui tmp = 0;
        ui v = C[i];
        ui tmp = sg.deg[deep][v] - sg.pIdx[v];
        // for(ui j = 0; j < edC; j++) {
        //     if(g.connectHash(C[i], C[j])) tmp++;
        // }
        // for(ui j = sg.pIdx[v]; j < sg.deg[deep][v]; j++) {
        //     ui w = sg.pEdge[j];
        //     if(C.isIn(w, 0, edC)) tmp++;
        // }

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }
#ifdef DDEBUG
if(uu == 2) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    if(pivotDeg+1 == edC && num == edC) {
#ifdef DDEBUG
d--;
#endif
        p += edC; updateAns(); return;
    }
// if(uu >= de) {
//     printf("p %u pd %u\n", pivot, pivotDeg); fflush(stdout);
// }
    C.changeTo(pivot, --edC);
    ui newEdC = 0;
    // for(ui i = 0; i < edC; i++) {
    //     if(g.connectHash(pivot, C[i])) 
    // }
    for(ui i = sg.pIdx[pivot]; i < sg.deg[deep][pivot]; i++) {
        ui w = sg.pEdge[i];
        if(C.isIn(w, 0, edC)) {
            C.changeTo(w, newEdC++);
        }
    }
    for(ui j = 0; j < newEdC; j++) {
        ui u = C[j];
        ui & ed = sg.deg[deep + 1][u];
        ed = sg.deg[deep][u];
        for(ui l = sg.pIdx[u]; l < ed; ) {
            ui w = sg.pEdge[l];
            if(!C.isIn(w, 0, newEdC)) {
                std::swap(sg.pEdge[l], sg.pEdge[--ed]);
            }
            else l++;
        }
    }

    searchSG(deep + 1, newEdC, p + 1, h);
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
        // for(ui j = 0; j < edC; j++) {
        //     if(g.connectHash(v, C[j])) C.changeToByPos(j, newEdC++); 
        // }
        for(ui j = sg.pIdx[v]; j < sg.deg[deep][v]; j++) {
            ui w = sg.pEdge[j];
            if(C.isIn(w, 0, edC)) {
                C.changeTo(w, newEdC++);
            }
        }
        for(ui j = 0; j < newEdC; j++) {
            ui u = C[j];
            ui & ed = sg.deg[deep + 1][u];
            ed = sg.deg[deep][u];
            for(ui l = sg.pIdx[u]; l < ed; ) {
                ui w = sg.pEdge[l];
                if(!C.isIn(w, 0, newEdC)) {
                    std::swap(sg.pEdge[l], sg.pEdge[--ed]);
                }
                else l++;
            }
        }

        searchSG(deep + 1, newEdC, p, h + 1);
    }

    freeMem(candSize);
#ifdef DDEBUG
d--;
#endif
}

void pivoter::search(ui edC, ui p, ui h) {
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
if(uu == 2) {
printf("pivot %u, pd %u\n", pivot, pivotDeg);
}
#endif
    if(pivotDeg+1 == edC && num == edC) {
#ifdef DDEBUG
d--;
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

    printf("pivoter.h\n");
}

