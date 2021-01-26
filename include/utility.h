#ifndef UTILITY_H
#define	UTILITY_H

#include <vector>
#include <deque>
#include <list>
#include <algorithm>

namespace Tailor
{
    template<class T> bool sign(T a, T b)
    {
        if (b >= 0) {
            return a;
        }
        else {
            return -a;
        }
    }

    template<class T> void remove_merge_dup(std::vector<T>& vec)
    {
        std::sort(vec.begin(), vec.end());

        auto beg = vec.begin();
        while(true)
        {
            beg = std::adjacent_find(beg, vec.end());
            auto temp_address = &(*beg);
            if (beg == vec.end()) return;

            beg->merge(*std::next(beg));
            vec.erase(std::next(beg));
            assert(&(*beg) == temp_address);
        }

        auto it = std::adjacent_find(vec.begin(), vec.end());
        assert(it == vec.end());
    }

    template<class T> void remove_merge_dup(std::deque<T>& vec)
    {
        std::sort(vec.begin(), vec.end());

        auto beg = vec.begin();
        while(true)
        {
            beg = std::adjacent_find(beg, vec.end());
            auto temp_address = &(*beg);
            if (beg == vec.end()) return;

            beg->merge(*std::next(beg));
            vec.erase(std::next(beg));
            assert(&(*beg) == temp_address);
        }

        auto it = std::adjacent_find(vec.begin(), vec.end());
        assert(it == vec.end());
    }

    template<class T> void remove_merge_dup(std::list<T>& con)
    {
        con.sort();

        auto beg = con.begin();
        while(true)
        {
            beg = std::adjacent_find(beg, con.end());
            auto temp_address = &(*beg);
            if (beg == con.end()) return;

            beg->merge(*std::next(beg));
            con.erase(std::next(beg));
            assert(&(*beg) == temp_address);
        }

        auto it = std::adjacent_find(con.begin(), con.end());
        assert(it == con.end());
    }

    template<class T> void remove_dup(std::vector<T>& con)
    {
        std::sort(con.begin(), con.end());
        con.erase(std::unique(con.begin(), con.end()), con.end());
    }

    template<class T> void remove_dup(std::deque<T>& con)
    {
        std::sort(con.begin(), con.end());
        con.erase(std::unique(con.begin(), con.end()), con.end());
    }

    template<class T> void remove_dup(std::list<T>& con)
    {
        con.sort();
        con.unique();
    }
}

#endif
