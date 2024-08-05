#include "kccPivot.h"
#include <algorithm>

// #define DDEBUG

double kccPivot::run() {
    printf("kccPivot.cpp\n");
    g.initHash();
    printf("initHash\n");

#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] > 0) {
        ui st = 0, ed = 0, edClique = 0;

        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], ed++);
        }
#ifdef DDEBUG
printf("node[0]:");
for(ui i = 0; i < ed; i++) printf("%u ", C[i]); printf("\n");
#endif

        ui maxDeg = 0, maxV = C[0];
        auto findM = [&](ui st, ui ed, ui & maxDeg, ui & maxV) {
            for(ui i = st; i < ed; i++) {
                ui deg = 0;
                for(ui j = st; j < ed; j++) {
                    if(g.connectHash(C[i], C[j])) deg++;
                }
                if(deg > maxDeg) {
                    maxDeg = deg;
                    maxV = C[i];
                }
            }
        };
        findM(st, ed, maxDeg, maxV);
        std::vector<ui> cliqueNei;
        edClique = ed;
        C.changeTo(maxV, --ed);
        ui tmpSt = 0;
        for(ui i = 0; i < ed; i++) {
            if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, tmpSt++);
        }

        while(tmpSt < ed) {
            ui maxDeg = 0, maxV = C[tmpSt];
            findM(tmpSt, ed, maxDeg, maxV);
            C.changeTo(maxV, --ed);

            ui newTmpSt = tmpSt;
            for(ui i = tmpSt; i < ed; i++) {
                if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, newTmpSt++);
            }
            tmpSt = newTmpSt;
        }

        if(edClique - ed >= k - 1) {
            answer += CN[edClique - ed][k - 1];
        }

        pClique[0].clear();
        for(ui i = ed; i < edClique; i++) pClique[0].push_back(C[i]);
        listing(0, ed, 0, 1);
    }

    return answer;
}

void kccPivot::listing(ui deep, ui edC, ui p, ui h) {
#ifdef DDEBUG
printf("    deep %u\n", deep);
printf("C:");
for(ui i = 0; i < edC; i++) printf("%u ", C[i]);printf("\n");
printf("clique:");
for(ui i = 0; i < pClique[deep].size(); i++) printf("%u ", pClique[deep][i]);printf("\n");
#endif

    if(h == k - 1) {
        // answer += p + edClique - st;
        // answer += p + edC + pClique[deep].size();
        answer += edC;
        return;
    }
    if(0 == edC) return;
    
    if(p + h + edC + pClique[deep].size() < k) return;
    
    ui pivot, pivotDeg = 0, pivotDegInClique;
    bool hasPivot = false;

    if(p == 0 || maxPreFixSize[deep] == 0) {
        for(ui i = 0; i < edC; i++) {
            ui deg = 0;
            for(ui j = 0; j < edC; j++) {
                if(g.connectHash(C[i], C[j])) deg++;
            }
            if(deg > pivotDeg) {
                hasPivot = true;
                pivotDeg = deg;
                pivot = C[i];

                ui degInClique = 0;
                for(auto v : pClique[deep]) {
                    if(g.connectHash(C[i], v)) degInClique++;
                }
                pivotDegInClique = degInClique;
            }
            else if(deg == pivotDeg){
                ui degInClique = 0;
                for(auto v : pClique[deep]) {
                    if(g.connectHash(C[i], v)) degInClique++;
                }
                if(degInClique > pivotDegInClique) {
                    pivotDegInClique = degInClique;
                    hasPivot = true;
                    pivotDeg = deg;
                    pivot = C[i];
                }
            }
        }

        ui newEdC = 0;
        for(ui i = 0; i < edC; i++) {
            if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++);
        }

        pClique[deep+1].clear();
        for(ui i = 0; i < pClique[deep].size(); i++) {
            if(g.connectHash(pivot, pClique[deep][i])) {
                pClique[deep+1].push_back(pClique[deep][i]);
            }
        }
        maxPreFixSize[deep+1] = std::max(maxPreFixSize[deep], (ui)pClique[deep+1].size());

        listing(deep + 1, newEdC, p+1, h);
    }
    else {
        Pset::iterator head = Po[deep].begin();
        Pset::iterator nextHead = Po[deep].begin();
        nextHead++;

        for(ui i = 0; i < edC; i++) {
            Pset::iterator l = head;
            Pset::iterator r = nextHead;
            ui u = C[i];
            ui cnt = 0;

            for(ui j = 0; j < head->first; j++) {
                if(!g.connectHash(u, pClique[deep][j])) {
                    cnt++; break;
                }
            }

            do {
                for(ui j = l->first; j < r->first; j++) {
                    if(!g.connectHash(u, pClique[deep][j])) {
                        cnt++; break;
                    }
                }
                l = r++;
            } while(cnt <= 1 && l != Po[deep].end());
            
            if(cnt > 1) continue;

            ui deg = 0;
            for(ui j = 0; j < edC; j++) {
                if(g.connectHash(C[i], C[j])) deg++;
            }
        }
    }

    if(hasPivot) {
        C.changeTo(pivot, st++);
        
        ui newSt = ed;
        for(ui i = st; i < newSt; ) {
            if(g.connectHash(pivot, C[i])) C.changeToByPos(i, --newSt);
            else i++;
        }

        ui candSize = newSt - st;
        ui * cand = allocMem(candSize);
        for(ui i = st; i < newSt; i++) cand[i-st] = C[i];

        listing(deep + 1, newSt, ed, edClique, p + 1, h);

        for(ui i = 0; i < candSize; i++) {
            C.changeTo(cand[i], st++);

            ui newSt = ed;
            for(ui j = st; j < newSt; ) {
                if(g.connectHash(cand[i], C[j])) C.changeToByPos(j, --newSt);
                else j++;
            }

            ui newEndClique = ed;
            for(ui j = ed; j < edClique; j++) {
                if(g.connectHash(cand[i], C[j])) C.changeToByPos(j, newEndClique++);
            }

            if(newEndClique-ed >= k-h-1) {
                answer += CN[newEndClique-ed][k-h-1];
            }
            listing(deep + 1, newSt, ed, newEndClique, p, h+1);
        }

        freeMem(candSize);

        return;
    }

    for(ui i = st; i < ed; i++) {
        ui u = C[i];
        C.changeTo(u, st++);

        ui newSt = ed;
        for(ui j = st; j < newSt; ) {
            if(g.connectHash(u, C[j])) C.changeToByPos(j, --newSt);
            else j++;
        }

        ui newEndClique = ed;
        for(ui j = ed; j < edClique; j++) {
            if(g.connectHash(u, C[j])) C.changeToByPos(j, newEndClique++);
        }

        if(newEndClique-ed >= k-h-1) {
            answer += CN[newEndClique-ed][k-h-1];
        }
        listing(deep + 1, newSt, ed, newEndClique, p, h+1);
    }
}

kccPivot::kccPivot(Graph && g, ui k) :g(g), k(k) {
    maxSize = g.coreNumber + 1;
    maxK = k + 1;
    answer = 0;
    computeC();

    initBuffer(g.coreNumber * (k+1));

    C.resize(g.n);

    pClique.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        pClique[i].resize(k);
    }

    preFixSize.resize(g.coreNumber);
    for(ui i = 0; i < g.coreNumber; i++) {
        preFixSize.resize(k);
    }
    maxPreFixSize.resize(g.coreNumber);

    Po.resize(g.coreNumber);

    printf("kccPivot2.cpp::kccPivot()\n");
}

void kccPivot::pivoter(ui deep, ui st, ui ed, ui p, ui h) {
    if(h == k - 1) {
        answer += p + ed - st;
        return;
    }
    if(ed == st) {
        answer += CN[p][k-h]; return;
    }
    if(ed == 1+st) {
        answer += CN[p+1][k-h]; return;
    }
    if(ed == 2+st) {
        if(g.connectHash(C[0], C[1])) answer += CN[p+2][k-h];
        else answer += CN[p+1][k-h] + CN[p][k-h-1];
        return;
    }

    ui pivot = C[st], pivotDeg = 0, num = 0;
    for(ui i = st; i < ed; i++) {
        ui tmp = 0;

        for(ui j = st; j < ed; j++) {
            if(g.connectHash(C[i], C[j])) tmp++;
        }

        if(tmp > pivotDeg) {
            pivot = C[i]; pivotDeg = tmp; num = 1;
        }
        else if(tmp == pivotDeg) num++;
    }
    if(pivotDeg+1+st == ed && num+st == ed) {
        answer += CN[p+ed-st][k-h]; return;
    }

    C.changeTo(pivot, --ed);
    ui newEdC = st;
    for(ui i = st; i < ed; i++) {
        if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++); 
    }
    pivoter(deep, st, newEdC, p+1, h);

    ui candSize = ed - st - pivotDeg;
    ui * cand = allocMem(candSize);
    memcpy(cand, C.begin() + pivotDeg, sizeof(ui) * candSize);

    for(ui i = 0; i < candSize; i++) {
        ui v = cand[i];
        C.changeTo(v, --ed);

        ui newEdC = st;
        for(ui j = st; j < ed; j++) {
            if(g.connectHash(v, C[j])) C.changeToByPos(j, newEdC++); 
        }

        pivoter(deep, st, newEdC, p, h + 1);
    }

    freeMem(candSize);
}