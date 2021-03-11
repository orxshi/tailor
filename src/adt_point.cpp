#include "adt.h"

namespace Tailor
{
    void ADTPoint::rotate(double angle, int axis, const Vector3& rot_point)
    {
        RotationMatrix rm;

        {
            Vector3 z;
            z(0) = dim_[0]-rot_point(0);
            z(1) = dim_[2]-rot_point(1);
            z(2) = dim_[4]-rot_point(2);
            auto newz = rm.rotate(angle, axis, z);

            dim_[0] = newz(0)+rot_point(0);
            dim_[2] = newz(1)+rot_point(1);
            dim_[4] = newz(2)+rot_point(2);
        }
        {
            Vector3 z;
            z(0) = dim_[1]-rot_point(0);
            z(1) = dim_[3]-rot_point(1);
            z(2) = dim_[5]-rot_point(2);
            auto newz = rm.rotate(angle, axis, z);

            dim_[1] = newz(0)+rot_point(0);
            dim_[3] = newz(1)+rot_point(1);
            dim_[5] = newz(2)+rot_point(2);
        }
    }

    void ADTPoint::move(const Vector3& v)
    {
        assert(dim_.size() == 4);

        dim_[0] += v(0);
        dim_[1] += v(0);
        dim_[2] += v(1);
        dim_[3] += v(1);
    }

    void ADTPoint::set_idx(int i)
    {
        idx_ = i;
    }

    ADTPoint::ADTPoint(): idx_(-1)
    {
    }

    ADTPoint::ADTPoint(const Vector3& p, int idx): idx_(idx)
    {
        dim_.resize(TAILOR_ADT_VAR);
        for (int i=0; i<TAILOR_ADT_DIM; ++i)
        {
            dim_[i*2]   = TAILOR_BIG_POS_NUM;
            dim_[i*2+1] = TAILOR_BIG_NEG_NUM;
        }

        for (int i=0; i<TAILOR_ADT_DIM; ++i)
        {
            dim_[i*2]   = std::min(dim_[i*2]  , p(i));
            dim_[i*2+1] = std::max(dim_[i*2+1], p(i));
        }
    }

    const std::vector<double>& ADTPoint::dim() const
    {
        return dim_;
    }

    double ADTPoint::dim(int i) const
    {
        assert(i >= 0);
        assert(i < dim_.size());

        return dim_[i];
    }

    int ADTPoint::idx() const
    {
        return idx_;
    }

    bool ADTPoint::overlap(const ADTPoint& other, bool verbose) const
    {
        //if (other.idx() == 537)
        //verbose = true;
        //else
        //verbose = false;

        for (int i=0; i<TAILOR_ADT_DIM; ++i)
        {
            if (verbose)
            {
                std::cout << "idx_ = " << idx_ << std::endl;

                std::cout << "dim_[0] = " << dim_[0] << std::endl;
                std::cout << "dim_[1] = " << dim_[1] << std::endl;
                std::cout << "dim_[2] = " << dim_[2] << std::endl;
                std::cout << "dim_[3] = " << dim_[3] << std::endl;
                std::cout << "dim_[4] = " << dim_[4] << std::endl;
                std::cout << "dim_[5] = " << dim_[5] << std::endl;

                std::cout << "other.dim(0) = " << other.dim(0) << std::endl;
                std::cout << "other.dim(1) = " << other.dim(1) << std::endl;
                std::cout << "other.dim(2) = " << other.dim(2) << std::endl;
                std::cout << "other.dim(3) = " << other.dim(3) << std::endl;
                std::cout << "other.dim(4) = " << other.dim(4) << std::endl;
                std::cout << "other.dim(5) = " << other.dim(5) << std::endl;
            }
            double dis1 = dim_[i*2] - other.dim(i*2+1);
            //double dis2 = dim_[i*2+1] - other.dim(i*2);
            double dis2 = other.dim(i*2) - dim_[i*2+1];
            if (verbose)
            {
                std::cout << "dis1 = " << dis1 << std::endl;
                std::cout << "dis2 = " << dis2 << std::endl;
                std::cout << "(dis1 > TAILOR_ZERO) = " << (dis1 > TAILOR_ZERO) << std::endl;
                std::cout << "(dis2 < TAILOR_ZERO) = " << (dis2 < TAILOR_ZERO) << std::endl;
            }
            //if (dim_[i*2] > other.dim(i*2+1)) return false;
            if (dis1 > TAILOR_ZERO) return false;
            //if (dis1 > 0.) return false;
            //if (dim_[i*2+1] < other.dim(i*2)) return false;
            //if (dis2 < TAILOR_ZERO) return false;
            if (dis2 > TAILOR_ZERO) return false;
            //if (dis2 < 0.) return false;
        }

        return true;
    }
}
