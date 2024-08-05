#include "kccPivot.h"
#include <algorithm>
#include <cassert>

// #define DDEBUG
// #define BASELINE

#define BRANCHES

#ifdef BRANCHES
ui HNODES = 0, PNODES = 0;
#endif

#ifdef DDEBUG
ui uu, uuu=2;
#endif

double kccPivot::run() {
    printf("kccPivot2.cpp\n");
    g.initHash();
    printf("initHash\n");

#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] > 0) {
        ui st = 0, ed = 0, edClique = 0;
#ifdef DDEBUG
uu = u;
#endif
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], ed++);
        }
#ifdef DDEBUG
if(uu==uuu){
printf("node[0]:");
for(ui i = 0; i < ed; i++) printf("%u ", C[i]); printf("\n");
}
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
#ifdef DDEBUG
double pans = answer;
#endif
        if(edClique - ed >= k - 1) {
            answer += CN[edClique - ed][k - 1];
        }

        pClique[0].clear();
        for(ui i = ed; i < edClique; i++) pClique[0].push_back(C[i]);
        listing(0, ed, 0, 1);
#ifdef DDEBUG
printf("perNode %u:%.0f\n", u, answer - pans);
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
    for(ui u = 0; u < g.n; u++) if((1<<u) & x) {
        for(ui v = u + 1; v < g.n; v++) if((1<<v) & x) {
            if(!g.connectHash(u, v)) return false;
        }
    }
    return true;
};

ui cnt = 0;
std::vector<ui> cnti(g.n);
for(uint32_t i = (1<<g.n)-1; i > 0; i--) {
    if(check(i)) {
        cnt++;
        for(ui u = 0; u < g.n; u++) if((1<<u) & i) {
            cnti[u]++; break;
        }
    }
}
printf("cnt:%u\n", cnt);
for(ui u = 0; u < g.n; u++) {
    printf("cnt%u:%u\n", u, cnti[u]);
}
#endif


#ifdef BRANCHES
printf("hnodes:%u\npNodes:%u\n", HNODES, PNODES);
#endif
    return answer;
}

void kccPivot::listing(ui deep, ui edC, ui p, ui h) {
#ifdef DDEBUG
if(uu==uuu){
printf("\n    deep %u, p%u, h%u, maxPsize %u\n", deep, p, h, maxPreFixSize[deep]);
printf("C:");
for(ui i = 0; i < edC; i++) printf("%u ", C[i]);printf("\n");
printf("clique:");
for(ui i = 0; i < pClique[deep].size(); i++) printf("%u ", pClique[deep][i]);printf("\n");
printf("Po:");
for(auto v:Po[deep]) printf("%u-%u ", v.first, v.second);printf("\n");
}
#endif
    auto updateAns = [&](ui a, ui b, ui t) {
#ifdef DDEBUG
if(uu==uuu){
printf("updateAns:%u %u %u, deep %u\n", a, b, t, deep);
}
#endif
        for(ui j = (t>b ? t-b : 0); j <= t && j <= a; j++) {
            answer += CN[a][j] * CN[b][t-j];
        }
    };

    if(h == k - 1) {
#ifdef DDEBUG
if(uu==uuu){
printf("updateAns edc:%u\n", edC);
}
#endif
        // answer += p + edClique - st;
        // answer += p + edC + pClique[deep].size();
        answer += edC;
        return;
    }
    if(0 == edC) return;
    
    if(p + h + edC + pClique[deep].size() < k) return;
    
    ui pivot, pivotDeg = 0, pivotDegInClique = 0, pivotIdx = 0;
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
            else if(deg == pivotDeg) {
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
#ifdef DDEBUG
if(uu==uuu){
if(hasPivot)
printf("pivot1 %u, pDeg %u, pDegInC %u\n", pivot, pivotDeg, pivotDegInClique);
}
#endif
    }
    else {
        for(ui i = 0; i < edC; i++) {
            ui idx = 0;
            Pset::iterator it = Po[deep].begin();
            ui u = C[i];
            ui cnt = 0;

            while(cnt < 1 && it != Po[deep].end()) {
                while(idx < it->first) {
                    if(!g.connectHash(u, pClique[deep][idx])) {
                        cnt++; 
                        idx = it->first;
                        break;
                    }
                    idx++;
                }
                ++it;
            }
            if(cnt == 1) {
                while(cnt == 1 && it != Po[deep].end()) {
                    while(idx < it->first) {
                        if(g.connectHash(u, pClique[deep][idx])) {
                            cnt++;
                            idx = it->first;
                            break;
                        }
                        idx++;
                    }
                    ++it;
                }
                if(cnt == 2) continue;

                while(idx < pClique[deep].size()) {
                    if(g.connectHash(u, pClique[deep][idx])) {
                        cnt++; break;
                    }
                    idx++;
                }
                if(cnt == 2) continue;
            }

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
#ifdef DDEBUG
if(uu==uuu){
if(hasPivot)
printf("pivot2 %u, pDeg %u, pDegInC %u\n", pivot, pivotDeg, pivotDegInClique);
}
#endif
    }


    auto hNode = [&](ui u) {
#ifdef BRANCHES
HNODES++;
#endif
        C.changeTo(u, --edC);

        std::vector<ui> & newClique = pClique[deep + 1];
        newClique.clear();
        Pset & nxtSet = Po[deep + 1];
        nxtSet.clear();

        ui i = 0, j = 1;
        Pset::iterator it = Po[deep].begin();

        while(it != Po[deep].end()) {
            while(i < it->first) {
                if(g.connectHash(u, pClique[deep][i])) {
                    newClique.push_back(pClique[deep][i]);
                }
                i++;
            }
            nxtSet.insert({newClique.size(), it->second});
            updateAns(p-j, newClique.size(), k-h-2);
            ++it;
            ++j;
        }
        maxPreFixSize[deep+1] = newClique.size();
        while(i < pClique[deep].size()) {
            if(g.connectHash(u, pClique[deep][i])) {
                newClique.push_back(pClique[deep][i]);
            }
            i++;
        }

        if(newClique.size() >= k-h-1) {
#ifdef DDEBUG
if(uu==uuu){
printf("updateAns clique+:C(%u, %u) %.0f\n", edC, CN[newClique.size()][k-h-1]);
}
#endif
            answer += CN[newClique.size()][k-h-1];
        }

        ui newEdC = 0;
        for(ui j = 0; j < edC; j++) {
            if(g.connectHash(u, C[j])) C.changeToByPos(j, newEdC++);
        }
        listing(deep + 1, newEdC, p, h+1);
    };

    if(hasPivot) {
        ui x = 0, i = 1;
        for(auto v: Po[deep]) {
            if(v.first < pivotDegInClique) {
                updateAns(p-i, v.first, k-h-2);
            }
            else {
                updateAns(p-i, pivotDegInClique, k-h-2);
            }
            i++;
        }
        updateAns(0, pivotDegInClique, k-h-1);

        C.changeTo(pivot, --edC);
        ui newEdC = 0;
        for(ui i = 0; i < edC; i++) {
            if(g.connectHash(pivot, C[i])) C.changeToByPos(i, newEdC++);
        }

        ui candSize = edC - newEdC;
        ui * cand = allocMem(candSize);
        for(ui i = newEdC; i < edC; i++) cand[i-newEdC] = C[i];

        pClique[deep + 1].clear();
        for(auto v : pClique[deep]) pClique[deep+1].push_back(v);
        Pset::iterator it = Po[deep].begin();
        i = 0;
        while(it != Po[deep].end()) {
            if(it->first >= pivotDegInClique) {
                ui w = it->first;
                for(ui j = i; j < w; ) {
                    if(g.connectHash(pivot, pClique[deep+1][j])) {
                        j++;
                    }
                    else std::swap(pClique[deep+1][j], pClique[deep+1][--w]);
                }
#ifdef DDEBUG
if(uu==uuu){
if(w != pivotDegInClique) {
    printf("w != pivotDegInClique: %u != %u\n pClique[d+1]:", w, pivotDegInClique);
    for(auto v: pClique[deep+1]) printf("%u ", v); printf("\n");
}
assert(w == pivotDegInClique);
}
#endif
                break;
            }
            i = it->first;
            ++it;
        }
        if(it == Po[deep].end()) {
            ui w = pClique[deep+1].size();
            for(ui j = i; j < w; ) {
                if(g.connectHash(pivot, pClique[deep+1][j])) {
                    j++;
                }
                else std::swap(pClique[deep+1][j], pClique[deep+1][--w]);
            }
        }
        // P.push_back(pivot);
        // preFixSize[]
        maxPreFixSize[deep+1] = std::max(maxPreFixSize[deep], pivotDegInClique);

        Po[deep + 1].clear();
        for(auto v: Po[deep]) Po[deep+1].insert(v);
        Po[deep + 1].insert({pivotDegInClique, pivot});
#ifdef BRANCHES
PNODES++;
#endif
        listing(deep + 1, newEdC, p + 1, h);


        for(ui i = 0; i < candSize; i++) {
            hNode(cand[i]);
        }

        freeMem(candSize);

        return;
    }
#ifdef DDEBUG
if(uu==uuu){
printf("noPivot deep %u\n", deep);
}
#endif
    while(edC) {
#ifdef DDEBUG
if(uu==uuu){
printf("deep %u, cand %u, pivot %u\n", deep, C[0], pivot);
}
#endif
        hNode(C[0]);
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

    // preFixSize.resize(g.coreNumber);
    // for(ui i = 0; i < g.coreNumber; i++) {
        // preFixSize.resize(k);
    // }
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