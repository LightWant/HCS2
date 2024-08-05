#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include <iostream>
#include <fstream>
#include <chrono>

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f_txt") && !ac.exist("-f") && !ac.exist("-f_txtD")
        && !ac.exist("-f_nd")) {
        return 0;
    }
    
    bool noVUM = false;
    if(ac.exist("noUVM")) noVUM = true;

    Graph g;
    if(ac.exist("-f_txt")) g.readFromText(ac["-f_txt"], noVUM);
    else if(ac.exist("-f_txtD")) g.readFromTextDoubleEdges(ac["-f_txtD"]);
    else if(ac.exist("-f_nd")) g.readFromTextNodes(ac["-f_nd"]);
    else g.readFromBin(ac["-f"]);


    std::cout << "load graph: n " << g.n << " m " << g.m << " maxD " << g.maxD << std::endl;

    auto s1 = std::chrono::steady_clock::now();

    g.changeToCoreOrder();

    auto s2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
    std::cout << "changeToCoreOrder:" << duration.count() << "ms" << std::endl;

    g.saveAsBin(ac["-d"]);

    // std::ofstream otdata(ac["-d"] + "data_noUVM.txt");
    // for(ui u = 0; u < g.n; u++) {
    //     for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
    //         otdata << u << ' ' << g.pEdge[i] << std::endl;
    //     }
    // }
    // otdata.close();

    std::ifstream in(ac["-label"]);
    std::ofstream ot(ac["-label"] + "_new");
    ui u, label;
    while(in >> u >> label) {
        ot << g.mp2[u] << ' ' << label << std::endl;
    }
    in.close();
    ot.close();

    return 0;
}