#include "cliqueLocal.h"

// #define DDEBUG
// #define BASELINE

void cliqueLocal::run() {
    // printf("cliqueLocal.cpp::run\n");
    g.initHash();
// g.print();
// for(ui i = 0; i < g.n; i++) {
//     printf("%u:", i);
//     for(ui j = pOIdx[i]; j < pOIdx[i+1]; j++) {
//         printf("%u ", pOEdge[j]);
//     }
//     printf("\n");
// }
#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) {
        ui edC = 0;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], edC++);
        }

        H[0] = u;
        search(edC, 0, 1);
    }

    // for(ui i = 0; i < g.n; i++) {
    //     printf("rid, now %u raw %u\n", i, g.mp[i]);
    // }
    // for(ui i = 0; i < g.m/2; i++) {
    //     ui u = 0;
    //     while(pOIdx[u + 1] <= reEegeId[i]) u++;
    //     printf("raw %u:%u-%u, now %u:%u-%u, %.0f\n", i, g.edges[i].first, g.edges[i].second, 
    //         reEegeId[i], u, pOEdge[reEegeId[i]],
    //         answers[reEegeId[i]]);
    // }

    for(ui i = 0; i < g.m/2; i++) {
        printf("%.0f\n", answers[reEegeId[i]]);
    }

#ifdef BASELINE
std::vector<ui> ec(g.m/2, 0);
auto print = [&](uint32_t x) {
    printf("c:");
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) printf("%u ", u);
    printf("\n");

    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        for(ui v = u + 1; v < g.n; v++) if((1<<v) & x) {
            auto st = pOEdge.begin() + pOIdx[u];
            auto ed = pOEdge.begin() + pOIdx[u + 1];
            ui idx = std::lower_bound(st, ed, v) - pOEdge.begin(); 
            ec[idx]++;
        }
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
for(ui e = 0; e < g.m/2; e++) printf("e_cnt %u, %u %.0f\n", e, ec[e], answers[e]);
fflush(stdout);
// for(ui e = 0; e < g.m/2; e++) assert(ec[e] == answers[e]);
#endif
}

#ifdef DDEBUG
ui deep = 0;
#endif
void cliqueLocal::search(ui edC, ui p, ui h) {
#ifdef DDEBUG
printf("\n");
printf("    deep %u, p %u, h %u\n", deep, p, h);
printf("C:");
for(ui i = 0; i < edC; i++) printf("%u ", C[i]); printf("\n");
printf("P:");
for(ui i = 0; i < p; i++) printf("%u ", P[i]); printf("\n");
printf("H:");
for(ui i = 0; i < h; i++) printf("%u ", H[i]); printf("\n");
printf("PP:"); for(auto i : PP) printf("%u ", i); printf("\n");
printf("PH:"); for(auto i : PH) printf("%u ", i); printf("\n");
printf("HH:"); for(auto i : HH) printf("%u ", i); printf("\n");
#endif

    auto updateAns = [&]() {
        if(k >= h+2 && p >= 2)
        for(auto i : PP) {
            answers[i] += CN[p-2][k-h-2];
        }
        if(k >= h+1 && p >= 1)
        for(auto i : PH) answers[i] += CN[p-1][k-h-1];
        for(auto i : HH) answers[i] += CN[p][k-h];
    };
    auto addEdges = [&](std::vector<ui> & edges, ui a, ui b) {
        if(a > b) std::swap(a, b);

        auto st = pOEdge.begin() + pOIdx[a];
        auto ed = pOEdge.begin() + pOIdx[a + 1];
        ui idx = std::lower_bound(st, ed, b) - pOEdge.begin(); 

        edges.push_back(idx);
    };
    
    if(h == k) {
        for(auto i : HH) answers[i] += 1;
        return;
    }
    if(edC == 0) {
        updateAns(); return;
    }
    if(edC == 1) {
        for(ui i = 0; i < p; i++) addEdges(PP, C[0], P[i]);
        for(ui i = 0; i < h; i++) addEdges(PH, C[0], H[i]);
        p += 1; 
    
        updateAns(); 

        p -= 1;
        for(ui i = 0; i < p; i++) PP.pop_back();
        for(ui i = 0; i < h; i++) PH.pop_back();
        return;
    }
    if(edC == 2) {
        if(g.connectHash(C[0], C[1])) {
            for(ui i = 0; i < p; i++) addEdges(PP, C[0], P[i]);
            for(ui i = 0; i < h; i++) addEdges(PH, C[0], H[i]);
            P[p] = C[0];
            p += 1; 
            for(ui i = 0; i < p; i++) addEdges(PP, C[1], P[i]);
            for(ui i = 0; i < h; i++) addEdges(PH, C[1], H[i]);
            p += 1; 
            
            updateAns();

            p -= 1; 
            for(ui i = 0; i < p; i++) PP.pop_back();
            for(ui i = 0; i < h; i++) PH.pop_back();
            p -= 1; 
            for(ui i = 0; i < p; i++) PP.pop_back();
            for(ui i = 0; i < h; i++) PH.pop_back();
        }
        else {
            for(ui i = 0; i < p; i++) addEdges(PP, C[0], P[i]);
            for(ui i = 0; i < h; i++) addEdges(PH, C[0], H[i]);
            P[p] = C[0];
            p += 1; 

            updateAns();

            p -= 1; 
            for(ui i = 0; i < p; i++) PP.pop_back();
            for(ui i = 0; i < h; i++) PH.pop_back();
//
            for(ui i = 0; i < p; i++) addEdges(PH, C[1], P[i]);
            for(ui i = 0; i < h; i++) addEdges(HH, C[1], H[i]);
            h += 1; 

            updateAns();

            h -= 1;
            for(ui i = 0; i < p; i++) PH.pop_back();
            for(ui i = 0; i < h; i++) HH.pop_back();
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

    if(pivotDeg+1 == edC && num == edC) {
        for(ui i = 0; i < edC; i++) {
            ui u = C[i];

            for(ui j = 0; j < p; j++) addEdges(PP, u, P[j]);
            for(ui j = 0; j < h; j++) addEdges(PH, u, H[j]);
            P[p++] = u;
        }

        updateAns(); 

        for(ui i = 0; i < edC; i++) {
            ui u = C[i];

            p--;
            for(ui j = 0; j < p; j++) PP.pop_back();
            for(ui j = 0; j < h; j++) PH.pop_back();
        }
        return;
    }
#ifdef DDEBUG
printf("deep %u, pivot %u\n", deep, pivot);
#endif

    C.changeTo(pivot, --edC);
    ui newEdC = 0;
    for(ui i = 0; i < edC; i++) {
        if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++); 
    }

    P[p] = pivot;
    for(ui i = 0; i < p; i++) addEdges(PP, pivot, P[i]);
    for(ui i = 0; i < h; i++) addEdges(PH, pivot, H[i]);
#ifdef DDEBUG
deep++;
#endif
    search(newEdC, p + 1, h);
#ifdef DDEBUG
deep--;
#endif
    for(ui i = 0; i < p; i++) PP.pop_back();
    for(ui i = 0; i < h; i++) PH.pop_back();

    ui candSize = edC - pivotDeg;
    ui * cand = allocMem(candSize);
    memcpy(cand, C.begin() + pivotDeg, sizeof(ui) * candSize);

    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];
        C.changeTo(v, --edC);
#ifdef DDEBUG
printf("deep %u, cand %u\n", deep, v);
#endif
        ui newEdC = 0;
        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(v, C[j])) C.changeToByPos(j, newEdC++); 
        }

        H[h] = v;
        for(ui j = 0; j < h; j++) addEdges(HH, v, H[j]);
        for(ui j = 0; j < p; j++) addEdges(PH, v, P[j]); 
#ifdef DDEBUG
deep++;
#endif
        search(newEdC, p, h + 1);
#ifdef DDEBUG
deep--;
#endif
        for(ui i = 0; i < h; i++) HH.pop_back();
        for(ui i = 0; i < p; i++) PH.pop_back();
    }

    freeMem(candSize);
}


cliqueLocal::cliqueLocal(Graph && g, ui k) :g(g), k(k) { 
    C.resize(g.n);

    maxSize = g.coreNumber + 10;
    maxK = k + 1;
    answers.resize(g.m / 2);

    pOEdge.resize(g.m / 2);
    pOIdx.resize(g.n + 1);
    for(ui u = 0; u < g.n; u++) {
        pOIdx[u + 1] = pOIdx[u];
        for(ui j = g.pIdx2[u]; j < g.pIdx[u + 1]; j++) {
            pOEdge[pOIdx[u+1]++] = g.pEdge[j];
        }
    }
    
    computeC();

    initBuffer(g.coreNumber * maxK);

    // P, H, PP, PH, HH
    P.resize(maxSize);
    H.resize(maxSize);

    ui maxSize2 = maxSize * maxSize;
    PP.resize(maxSize2);
    PH.resize(maxSize2);
    HH.resize(maxSize2);
    PP.clear();
    PH.clear();
    HH.clear();

//only when g.edges is sorted
    reEegeId.resize(g.m/2);
    for(ui i = 0; i < g.m/2; i++) {
        ui u = g.edges[i].first;
        ui v = g.edges[i].second;
        u = g.mp2[u];
        v = g.mp2[v];
        if(u > v) std::swap(u, v);

        auto st = pOEdge.begin() + pOIdx[u];
        auto ed = pOEdge.begin() + pOIdx[u + 1];
        ui idx = std::lower_bound(st, ed, v) - pOEdge.begin(); 

        reEegeId[i] = idx;
    }

    // printf("cliqueLocal.h\n");
}
