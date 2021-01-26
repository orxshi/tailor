#ifndef ADT_HPP
#define	ADT_HPP

namespace Tailor
{
    template<typename Rai> ADTPoint::ADTPoint(Rai begin, Rai end, int idx): idx_(idx)
    {
        dim_.resize(TAILOR_ADT_VAR);
        for (int i=0; i<TAILOR_ADT_DIM; ++i)
        {
            dim_[i*2]   = TAILOR_BIG_POS_NUM;
            dim_[i*2+1] = TAILOR_BIG_NEG_NUM;
        }

        for (auto it=begin; it!=end; ++it)
        {
            for (int i=0; i<TAILOR_ADT_DIM; ++i)
            {
                dim_[i*2]   = std::min(dim_[i*2]  , it->r(i));
                dim_[i*2+1] = std::max(dim_[i*2+1], it->r(i));
            }
        }
    }
}

#endif
