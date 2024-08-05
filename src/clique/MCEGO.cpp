#include "kccPivot.h"

ui kccPivot::maxiDeg(ui st, ui ed, ui xx) {
    //find from x vertices with maximum degree
    constexpr ui x = 5;
// assert(ed - st > x);
    ui maxDegX[x+1] = {0}, maxVX[x+1] = {0};
    ui j = 0;
    for(ui i = st; i < st+xx; i++) {
        ui deg = 0;
        for(ui j = st; j < ed; j++) {
            if(g.connectHash(C[i], C[j])) deg++;
        }
        maxVX[i-st+1] = C[i];
        maxDegX[i-st+1] = deg;
    }
    for(ui i = st+xx; i < ed; i++) {
        ui deg = 0;
        for(ui j = st; j < ed; j++) {
            if(g.connectHash(C[i], C[j])) deg++;
        }
        if(deg > maxDegX[xx]) {
            ui j = xx-1;
            for(; j >= 1; j--) if(maxDegX[j] >= deg) break;
            j++;
            do {
                std::swap(maxDegX[j++], deg);
            } while(j <= xx);
        }
    }

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

    auto solve = [&](ui maxV, ui maxDeg) {
        ui stClique = ed;
        C.changeTo(maxV, --stClique);
        ui tmpSt = 0;
        for(ui i = 0; i < stClique; i++) {
            if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, tmpSt++);
        }

        while(tmpSt < stClique) {
            ui maxDeg = 0, maxV = C[tmpSt];
            findM(tmpSt, ed, maxDeg, maxV);
            C.changeTo(maxV, --stClique);

            ui newTmpSt = tmpSt;
            for(ui i = tmpSt; i < stClique; i++) {
                if(!g.connectHash(maxV, C[i])) C.changeToByPos(i, newTmpSt++);
            }
            tmpSt = newTmpSt;
        }

        return stClique;
    };

    ui minStClique = solve(maxVX[1], maxDegX[1]);
    ui tt = 1;
    for(ui t = 2; t <= xx; t++) {
        ui stClique = solve(maxVX[t], maxDegX[t]);
        if(stClique < minStClique) {
            minStClique = stClique;
            tt = t;
        }
    }

    return solve(maxVX[tt], maxDegX[tt]);
}

ui kccPivot::coreClique(ui st, ui ed, ui core) {
    //core reduction


    return 0;
}

ui kccPivot::MCEGO(ui st, ui ed) {
    //maximum degree-based heuristic
    ui stClique = maxiDeg(st, ed);
    ui cliqueSize = ed - stClique;

    st = coreClique(st, ed, cliqueSize);

    return coreClique(st, ed, cliqueSize);
}