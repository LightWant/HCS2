#include "kccPivot.h"
#include <algorithm>
#include <cassert>

// #define DDEBUG

#define BRANCHES

#ifdef BRANCHES
ui HNODES = 0, PNODES = 0;
#endif

double kccPivot::run() {
    printf("kccPivot.cpp\n");
    g.initHash();
    printf("initHash\n");

#ifdef DDEBUG
g.print();
#endif

    for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] > 0) {
        ui st = 0, ed = 0, edClique = 0;
#ifdef DDEBUG
printf("    st %u\n", u);
#endif
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            C.changeTo(g.pEdge[i], ed++);
        }
#ifdef DDEBUG
printf("node[0]:");
for(ui i = 0; i < ed; i++) printf("%u ", C[i]); printf("\n");
#endif

        // ui maxDeg = 0, maxV = C[0];
        // auto findM = [&](ui st, ui ed, ui & maxDeg, ui & maxV) {
        //     for(ui i = st; i < ed; i++) {
        //         ui deg = 0;
        //         for(ui j = st; j < ed; j++) {
        //             if(g.connectHash(C[i], C[j])) deg++;
        //         }
        //         if(deg > maxDeg) {
        //             maxDeg = deg;
        //             maxV = C[i];
        //         }
        //     }
        // };
        // findM(st, ed, maxDeg, maxV);
        // std::vector<ui> cliqueNei;
        // edClique = ed;
        // C.changeTo(maxV, --ed);
        // ui tmpSt = 0;
        // for(ui i = 0; i < ed; i++) {
        //     if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, tmpSt++);
        // }

        // while(tmpSt < ed) {
        //     ui maxDeg = 0, maxV = C[tmpSt];
        //     findM(tmpSt, ed, maxDeg, maxV);
        //     C.changeTo(maxV, --ed);

        //     ui newTmpSt = tmpSt;
        //     for(ui i = tmpSt; i < ed; i++) {
        //         if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, newTmpSt++);
        //     }
        //     tmpSt = newTmpSt;
        // }

        edClique = ed;

        if(ed - st > 20) ed = maxiDeg(st, ed, 5);
        else ed = maxiDeg(st, ed, 1);
// for(ui i = edClique; i < ed; i++) {
//     for(ui j = i + 1; j < ed; j++) {
//         assert(g.connectHash(C[i], C[j]));
//     }
// }
        if(edClique - ed >= k - 1) {
            answer += CN[edClique - ed][k - 1];
        }
#ifdef DEBUG
double pans = answer, ans = 0;
for(ui i = st; i < ed; i++) {
    for(ui j = i + 1; j < edClique; j++) {
        if(g.connectHash(C[i], C[j])) ans++;
    }
}
#endif
        listing(0, st, ed, edClique, 0, 1);
#ifdef DEBUG
if(pans + ans != answer) {
    printf("error %u: %.0f %.0f %.0f\n", u, pans, ans, answer);
}
#endif
    }
#ifdef BRANCHES
printf("hnodes:%u\npNodes:%u\n", HNODES, PNODES);
#endif
    return answer;
}

void kccPivot::listing(ui deep, ui st, ui ed, ui edClique, ui p, ui h) {
#ifdef DDEBUG
printf("    deep %u, p %u, h %u\n", deep, p, h);
printf("C:");
for(ui i = st; i < ed; i++) printf("%u ", C[i]);printf("\n");
printf("clique:");
for(ui i = ed; i < edClique; i++) printf("%u ", C[i]);printf("\n");

#endif

    if(h == k-1) {
        // answer += p + edClique - st;
        answer += ed - st;
// #ifdef DDEBUG
// printf("ans1 += %u+%u\n", p, ed-st);
// #endif
        return;
    }
    if(st == ed) {   
//         if(p >= k-h) {
// #ifdef DDEBUG
// printf("ans2 += C(%u, %u)\n", p, k-h);
// #endif
//             answer += CN[p][k - h];
//         }
        return;
    }

    // if(st+1 == ed) {
    //     ui neiInClique = 0;
    //     for(ui j = ed; j < edClique; j++) {
    //         if(g.connectHash(C[st], C[j])) {
    //             neiInClique++;
    //         }
    //     }
    //     if(neiInClique == 0) {
    //         if(p+1 >= k-h) {
    //             answer += CN[p+1][k-h];
    //         }
    //     }
    //     else {
    //         if(neiInClique >= k-h-1) answer += CN[neiInClique][k-h-1];
    //         if(p >= k-h-1) answer += CN[p][k-h-1];
    //     }
        
    //     return;
    // }


    // if(p + h + edClique - st < k) return;

    // if(edClique == ed) {
    //     pivoter(deep, st, ed, p, h);
    //     return;
    // }
    
    ui pivot = C[st], pivotDeg = 0;
    bool hasPivot = false;
    for(ui i = st; i < ed; i++) {
        bool nonNeiInClique = true;
        for(ui j = ed; j < edClique; j++) {
            if(g.connectHash(C[i], C[j])) {
                nonNeiInClique = false;
                break;
            }
        }
        if(!nonNeiInClique) continue;

        ui deg = 0;
        for(ui j = st; j < ed; j++) {
            if(g.connectHash(C[i], C[j])) deg++;
        }

        if(!hasPivot) {
            hasPivot = true;
            pivot = C[i];
            pivotDeg = deg;
        }
        else if(deg > pivotDeg) {
            pivotDeg = deg;
            pivot = C[i];
        }
    }
#ifdef DDEBUG
if(hasPivot) {
    printf("hasPivot!\n");
    printf("pivot %u\n", pivot);
}
else {
    printf("do not have pivot\n");
}
#endif
    // if(hasPivot) hasPivot = false;
// if(st+1==ed) assert(hasPivot == false);

// if(st+1==ed && hasPivot) return;

    if(hasPivot) {
// for(ui i = ed; i < edClique; i++) assert(!g.connectHash(pivot, C[i]));
        C.changeTo(pivot, st++);
        
        ui newSt = ed;
        for(ui i = st; i < newSt; ) {
            if(g.connectHash(pivot, C[i])) C.changeToByPos(i, --newSt);
            else i++;
        }
// assert(ed > newSt);
        
// for(ui i = newSt; i < ed; i++) assert(g.connectHash(pivot, C[i]));
// for(ui i = st; i < newSt; i++) assert(!g.connectHash(pivot, C[i]));
// assert(C[st-1] == pivot);
        if(p >= k-h-1) {
            answer += CN[p][k-h-1];
        }
#ifdef BRANCHES
PNODES++;
#endif
        listing(deep + 1, newSt, ed, edClique, p + 1, h);

        ui candSize = newSt - st;
        if(candSize == 0) return;
        ui * cand = allocMem(candSize);
        for(ui i = st; i < newSt; i++) cand[i-st] = C[i];
#ifdef BRANCHES
HNODES += candSize;
#endif
        for(ui i = 0; i < candSize; i++) {
            C.changeTo(cand[i], st++);
// assert(!g.connectHash(pivot, cand[i]));
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
            if(p >= k-h-1) {
                answer += CN[p][k-h-1];
            }
#ifdef DDEBUG
printf("deep %u, cand_non_pivot %u\n", deep, cand[i]);
#endif
            listing(deep + 1, newSt, ed, newEndClique, p, h+1);
        }

        freeMem(candSize);

        return;
    }
#ifdef BRANCHES
HNODES += ed - st;
#endif
    for(ui i = st; i < ed; i++) {
        ui u = C[i];
        st++;
        // C.changeTo(u, st++);
#ifdef DDEBUG
printf("deep %u, cand_non_pivot %u\n", deep, u);
#endif
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
        if(p >= k-h-1) {
            answer += CN[p][k-h-1];
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

    printf("kccPivot.cpp::kccPivot()\n");
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