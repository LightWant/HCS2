#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

void printUsage() {
    std::cout << "bin/run " << std::endl;
    std::cout << "-f graph file directory(edge.bin & idx.bin)" << std::endl;
}

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f")) {
        printUsage();
        return 0;
    }

    Graph g;
    g.readFromBin(ac["-f"]);

    //std::cout << "load graph: n " << g.n << " m " << g.m << " maxD " << g.maxD << std::endl;

    for(ui u = 0; u < g.n; u++) {
        for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            printf("%u %u\n", u, g.pEdge[i]);
        }
    }

    return 0;
}