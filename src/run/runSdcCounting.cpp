#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../sdc/sdcCounting2DP.h"
#include <iostream>
#include <chrono>
#include <iomanip>

void printUsage() {
    std::cout << "bin/run " << std::endl;
    std::cout << "-f graph file directory(edge.bin & idx.bin)" << std::endl;
    std::cout << "-f_txt graph file text file, each edge exists one time" << std::endl;
    std::cout << "-f_txtD graph file text file, each edge exists two times" << std::endl;
    std::cout << "-s" << std::endl;
    std::cout << "-q" << std::endl;
}

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f_txt") && !ac.exist("-f") && !ac.exist("-f_txtD")) {
        printUsage();
        return 0;
    }
    else if(!ac.exist("-s") || !ac.exist("-q")) {
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
// for(ui i = g.pIdx[0]; i < g.pIdx[1]; i++) {
//     printf("%u ", g.pEdge[i]);
// }printf("\n");fflush(stdout);
    g.changeToCoreOrder();
    for(ui i = 0; i < g.n; i++) assert(g.pIdx[i+1] - g.pIdx2[i] <= g.coreNumber);

    auto s2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
    std::cout << "changeToCoreOrder:" << duration.count() << "ms" << std::endl;
    std::cout << "coreNumber:" << g.coreNumber << std::endl;

    int s = std::atoi(ac["-s"].c_str());
    int q = std::atoi(ac["-q"].c_str());
    int Q = 20;
    if(ac.exist("-Q")) Q = std::stoi(ac["-Q"]);
    std::cout << "s:" << s << " q:" << q << " Q:" << Q << std::endl;

    sdcCounting pS(std::move(g), s, q, Q);

    auto t1 = std::chrono::steady_clock::now();
    
    auto cnt = pS.run();

    auto t2 = std::chrono::steady_clock::now();
    auto durationt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << durationt.count() << "ms" << std::endl;

    for(ui i = q; i < cnt.size(); i++) {
        if(cnt[i] <= 0.5) break;
        std::cout << std::fixed << std::setprecision(0) <<
             i << "-sdc:" << cnt[i] << std::endl;
    }
    return 0;
}