#ifndef CANTOR_H
#define CANTOR_H

#include <utility>
#include <cmath>

namespace Tailor
{
    int pair(int x, int y);
    std::pair<int, int> unpair(int z);
    //int pair(int x, int y, int t);
    //void unpair(int z, int& x, int& y, int& t);

    int cantor(int x, int y);
    std::pair<int, int> unpair_cantor(int z);
    int szudzik(int x, int y);
    std::pair<int, int> unpair_szudzik(int z);
}

#endif
