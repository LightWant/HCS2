#include "candidate.hpp"
#include <iostream>

int main()
{
    unsigned n = 10, k = 5;
    candidate c(n, k);

    // c.insert(6, 3);
    // c.insert(4, 3);
    // c.insert(7, 4);
    // c.insert(2, 3);
    c.insert(1, 0);
    c.insert(7, 0);

    std::cout << "inserted" << std::endl;

    for(auto v = c.ibegin(); v != c.iend(); v = c.toNxt(v)) {
        std::cout << v << std::endl;
    }

    c.remove(6);
    c.remove(7);
    c.insert(6, 2);

    std::cout << "modified" << std::endl;

    for(auto v = c.ibegin(); v != c.iend(); v = c.toNxt(v)) {
        std::cout << v << std::endl;
    }


    std::cout << "ranged for" << std::endl;
    for(auto v : c) std::cout << v << std::endl;

    return 0;
}