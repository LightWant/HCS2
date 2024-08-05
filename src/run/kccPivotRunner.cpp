#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../clique/kccPivot.h"
#include <iostream>
#include <chrono>
#include <iomanip>

void printUsage() {
    std::cout << "bin/run " << std::endl;
    std::cout << "-f graph file directory(edge.bin & idx.bin)" << std::endl;
    std::cout << "-f_txt graph file text file, each edge exists one time" << std::endl;
    std::cout << "-f_txtD graph file text file, each edge exists two times" << std::endl;
    std::cout << "-k" << std::endl;
}

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f_txt") && !ac.exist("-f") && !ac.exist("-f_txtD")) {
        printUsage();
        return 0;
    }
    
    bool noVUM = false;
    if(ac.exist("noUVM")) noVUM = true;

    Graph g;
    if(ac.exist("-f_txt")) g.readFromText(ac["-f_txt"], noVUM);
    else if(ac.exist("-f_txtD")) g.readFromTextDoubleEdges(ac["-f_txtD"]);
    else g.readFromBin(ac["-f"]);


    std::cout << "load graph: n " << g.n << " m " << g.m << " maxD " << g.maxD << std::endl;

    auto s1 = std::chrono::steady_clock::now();
// ui maxD = 0;
// for(ui u = 0; u < g.n; u++) {
//     maxD = std::max(maxD, g.pIdx[u + 1] - g.pIdx2[u]);
// }
// std::cout << "coreNumber:" << maxD << std::endl;
    g.changeToCoreOrder();
    // for(ui i = 0; i < g.n; i++) assert(g.pIdx[i+1] - g.pIdx2[i] <= g.coreNumber);

    auto s2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
    std::cout << "changeToCoreOrder:" << duration.count() << "ms" << std::endl;
    std::cout << "coreNumber:" << g.coreNumber << std::endl;

    int k = 7;
    if(ac.exist("-k")) k =  std::stoi(ac["-k"]);
    std::cout << "k:" << k << std::endl;

    kccPivot pP(std::move(g), k);

    auto t1 = std::chrono::steady_clock::now();
    
    auto cnt = pP.run();

    auto t2 = std::chrono::steady_clock::now();
    auto durationt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << durationt.count() << "ms" << std::endl;

    std::cout << std::fixed << std::setprecision(0) <<
             k << "-clique:" << cnt << std::endl;
    
    return 0;
}