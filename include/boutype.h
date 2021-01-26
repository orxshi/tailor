#ifndef TAILOR_BOUTYPE_H
#define TAILOR_BOUTYPE_H

namespace Tailor
{
    enum class BouType
    {
        undefined = -1,
        wall = 1,
        dirichlet = 2,
        empty = 3,
        interior = 4,
        farfield = 9,
        partition = 10,
        interog = 11,
    };
}

#endif
