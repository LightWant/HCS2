#include "d2kraw.h"

ui d2kRaw::run() {
    std::vector<bool> vis(g.n, false);
    vc nodesInCX;
    vc cnn(g.n, 0);

    for(ui u = 0; u < g.n; u++) {
#ifdef DDDEBUG
std::cout<<"    start "<<u<<' '<<answer<<std::endl;
#endif
        //P is empty
        //X is 2-hop neighbors < u
        //C is 2-hop neighbors > u
        if(g.pIdx[u + 1] - g.pIdx2[u] + k < q) continue;
    
        vc C;
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) 
            C.push_back(g.pEdge[i]);
        
        ui sz = 0;
        do {
            sz = C.size();
            vc filtered;
            for(ui i = 0; i < sz; i++) {
                ui d = 0;
                for(ui j = 0; j < sz; j++) {
                    if(g.connect(C[i], C[j])) d++;
                }
                if(d + 2*k >= q) filtered.push_back(C[i]);
            }
            C = filtered;
        } while(C.size() != sz);

        vc X;
        for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++)
            X.push_back(g.pEdge[i]);
        do {
            sz = X.size();
            vc filtered;
            for(ui i = 0; i < sz; i++) {
                ui d = 0;
                for(ui j = 0; j < C.size(); j++) {
                    if(g.connect2(X[i], C[j])) d++;
                }
                if(d + 2*k >= q) filtered.push_back(X[i]);
            }
            X = filtered;
        } while(X.size() != sz);

        for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) vis[g.pEdge[i]] = true;
        vis[u] = true;
        sz = C.size();
        for(ui i = 0; i < sz; i++) {
            ui v = C[i];
            for(ui j = g.pIdx[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if(vis[w]) continue;

                ui d = 0;
                for(ui l = 0; l < sz; l++) if(g.connect2(w, C[l])) d++;

                if(d + 2*k >= q+2) {
                    if(w > u) C.push_back(w);
                    else X.push_back(w);
                    vis[w] = true;
                    addNN(cnn, w, u);
                    nodesInCX.push_back(w);
                }
            }
        }
        for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) vis[g.pEdge[i]] = false;
        vis[u] = false;
        for(auto v :  nodesInCX) vis[v] = false;
        nodesInCX.clear();
        
        vc P(C.size()); P.clear();
        P.push_back(u);

        //build sub-graph g
        auto buildV = [&](ui v) {
            sg.pIdx[v] = sg.pIdx2[v] = g.pIdx[v];
#ifdef DDDEBUG
printf("%u:", v);
#endif
            for(ui i = 0; i < C.size(); i++) {
                if(g.connect(v, C[i])) {
                    sg.pEdge[sg.pIdx2[v]++] = C[i];
#ifdef DDDEBUG
printf("%u ", C[i]);
#endif
                }
            }
#ifdef DDDEBUG
printf("\n");
#endif
        };

        for(ui i = 0; i < C.size(); i++) buildV(C[i]);
        for(ui i = 0; i < X.size(); i++) buildV(X[i]);
        sg.pIdx[u] = sg.pIdx2[u] = g.pIdx[u];
        for(ui i = 0; i < sz; i++) sg.pEdge[sg.pIdx2[u]++] = C[i];
        
        bkPivot(0, P, C, X, cnn);
    }

    return answer;
}

void d2kRaw::bkPivot(ui deep, vc & P, vc & C, vc & X, vc & cnn) {
    if(C.size() == 0) {
#ifdef DDDEBUG
printf("    edC == 0\n");
#endif
        if(X.size() == 0 && P.size() >= q) {
#ifdef DDDEBUG
printf("    P is ans!\n");
#endif
            answer++;
        }
        return;
    }

    if(C.size() + P.size() < q) return;

    //is P+C a k-plex
    ui minDPC = g.n;
    //find pivot in C and X with max degree in C
    ui maxD = 0, pv = C[0];
    bool inC = true;


}
