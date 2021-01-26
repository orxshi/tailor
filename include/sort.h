#ifndef SORT_H
#define	SORT_H

#include <algorithm>

namespace Tailor
{
    template<class Iter>
        void merge_sort(Iter first, Iter last)
        {
            if (last - first > 1)
            {
                Iter middle = first + (last - first) / 2;
                merge_sort(first, middle); // [first, middle)
                merge_sort(middle, last);  // [middle, last)
                std::inplace_merge(first, middle, last);
            }
        }
}

#endif
