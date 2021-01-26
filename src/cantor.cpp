#include "cantor.h"

namespace Tailor
{
    //int pair(int x, int y, int t)
    //{
    //    int z = cantor(x, y);
    //    return cantor(z, t);
    //}

    //void unpair(int z, int& x, int& y, int& t)
    //{
    //    auto m = unpair_cantor(z);
    //    t = m.second;
    //    auto n = unpair_cantor(m.first);
    //    x = n.first;
    //    y = n.second;
    //}

    int pair(int x, int y)
    {
        return cantor(x, y);
    }

    std::pair<int, int> unpair(int z)
    {
        return unpair_cantor(z);
    }

    int cantor(int x, int y)
    {
        int first, second;
        if (x < y)
        {
            first = x;
            second = y; 
        }
        else
        {
            first = y;
            second = x;
        }

        int z =  (x + y + 1) * (x + y) / 2 + y;

        return z;
    }

    std::pair<int, int> unpair_cantor(int z)
    {
        //https://math.stackexchange.com/a/222835/142155

        int n = 0.5 * (std::sqrt(1 + 8 * z) - 1);
        int tn = 0.5 * n * (n + 1);
        int y = z - tn;
        int x = n - y;

        return std::pair<int, int>(x, y);
    }

    int szudzik(int x, int y)
    {
        // https://www.vertexfragment.com/ramblings/cantor-szudzik-pairing-functions/

        int z;
        if (x < y) {
            z = std::pow(y, 2) + x;
        }
        else {
            z = std::pow(x, 2) + x + y;
        }

        return z;
    }

    std::pair<int, int> unpair_szudzik(int z)
    {
        // https://gist.github.com/antimatter15/8cb2538f4bd195e0b439560ec8c8e5b9

        int q = std::floor(std::sqrt(z));
        int l = z - std::pow(q, 2);
        if (l < q) {
            return std::pair<int, int>(l, q);
        }
        else {
            return std::pair<int, int>(q, l-q);
        }
    }
}
