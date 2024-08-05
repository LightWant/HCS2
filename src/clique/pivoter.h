#ifndef PIVOTER_H
#define PIVOTER_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <vector>
#include <algorithm>
#include <set>
#include <cassert>

class pivoter
{
private:
    Graph g;
    ui k = 0;
    LinearSet C;
    std::vector<double> answers;
    
private://combinatorial number
    ui maxK = 31;
    ui maxSize = 3000;
    double ** CN = nullptr, *bf3 = nullptr;
    void computeC() {
        ui maxD = maxSize, maxD2 = maxK;
        CN = new double*[maxD];
        bf3 = new double[maxD2 * maxD];
        for(int i = 0; i < maxD; i++) {
            CN[i] = bf3 + i * maxD2;
        }
        CN[0][0] = 1;
        CN[1][0] = 1;
        CN[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            CN[i][0] = 1;
            if(i < maxD2) CN[i][i] = 1;
            for(int j = 1; j < i && j < maxD2; j++) {
                CN[i][j] = CN[i - 1][j - 1] + CN[i - 1][j];
            }
        }
    }

private://raw pivoter
    void search(ui edC, ui p, ui h);
    
private://sg + kclist
    Graph sg;
    void searchSG(ui deep, ui edC, ui p, ui h);
    //real sg(pivoterSG.cpp)
    void searchSG(ui deep, ui p, ui h);

private://kclist
    std::vector<std::vector<ui>> nodes;
    std::vector<ui> level;
    std::vector<bool> ok;
    void searchKCList(ui deep, ui p, ui h);

private://pivoter + find a large clique + kclist
//st-ed, stC-edC: cands, stC-edC:a clique
    void searchSGClique(ui st, ui ed, ui stC, ui edC, ui p, ui h);
    ui expandClique(ui stC, ui edC);
    void findTheMax(ui st, ui ed, ui & maxDeg, ui & maxV);

private://small pivot
    void search(ui stC, ui edC, ui p, ui h);

private:
//sgClique, extract a large clique at first, then build the 
//search tree.
    std::vector<ui> clique, Ps, Hs;
    LinearSet candL, candR;
    std::vector<std::vector<ui>> tmpNodesL, tmpNodesR;
    void searchSGClique(ui deep, ui edC, ui p, ui h);
    struct treePath {
        int p1, h1, p2, h2;
    };
    std::vector<double> ansBiclique; 
    //ansBiclique[i] is the count of (p,q)-biclique that p+q=i
    ui limitOfPQ = 0;
    void pivotCount(int l, int pL, int pR, treePath t);
    struct c4{
        constexpr static ui sz = 4;
        ui v[sz];
        // std::vector<ui> v;
        c4() {}
        c4(ui a,ui b, ui c, ui d) {
            v[0] = a;
            v[1] = b;
            v[2] = c;
            v[3] = d;
            std::sort(v, v + sz);
        }
        bool isClique(Graph * g) {
            for(ui i = 0; i < sz; i++) {
                for(ui j = i + 1; j < sz; j++) {
                    if(!g->connectHash(v[i], v[j])) return false;
                }
            }
            return true;
        }
        void sort() {
            std::sort(v, v + 4);
        }
        bool operator ==(const c4 & t) const {
            for(ui i = 0; i < sz; i++) {
                if(v[i] != t.v[i]) return false;
            }
            return true;
        }
        bool operator < (const c4 & t) const {
            if(*this == t) return false;
            for(ui i = 0; i < sz; i++) {
                if(v[i] != t.v[i]) return v[i] < t.v[i];
            }
            return true;
        }
    };
    std::set<c4> sc4;
    

private:
    ui * memBuffer = nullptr;
    ui * pointerBuffer = nullptr;
    void initBuffer(ui n) {
        memBuffer = new ui[n];
        pointerBuffer = memBuffer;
    }
    ui * allocMem(ui n) {
        ui * tmp = pointerBuffer;
        pointerBuffer += n;
        return tmp;
    }
    void freeMem(ui n) {
        pointerBuffer -= n;
    }

public:
    pivoter(Graph && g, ui k);
    ~pivoter() { 
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
        if(memBuffer != nullptr) delete [] memBuffer; 
    }

    std::vector<double> run();
};


#endif 