#ifndef ROTATIONMATRIX_H
#define	ROTATIONMATRIX_H

//#include <vec3.h>
#include "matrix.h"

namespace Tailor
{
    class RotationMatrix
    {
        public: 

        Vector3 rotate(double angle, int axis, const Vector3& v);

        private:

        using trimat = std::array<double, 9>;
        trimat rx, ry, rz;

        void set_rx(double yaw);
        void set_ry(double pitch);
        void set_rz(double roll);
        Vector3 mul(const trimat& m, const Vector3& v);
    };
}

#endif
