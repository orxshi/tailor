#include "matrix.h"

namespace Tailor
{
    Vector5 operator/(const Vector5& f, const Matrix5& A)
    {
        // perform x = inv (A) * f

        double b11,b21,b22,b31,b32,b33,b41,b42,b43,b44,b51,b52,b53,b54,b55;
        double u12,u13,u14,u15,u23,u24,u25,u34,u35,u45;
        double d1,d2,d3,d4,d5;

        // decompose A into L and U
        b11 = 1 / A(0,0);
        u12 = A(0,1) * b11;
        u13 = A(0,2) * b11;
        u14 = A(0,3) * b11;
        u15 = A(0,4) * b11;
        b21 = A(1,0);
        b22 = 1 / (A(1,1) - b21 * u12);
        u23 = (A(1,2) - b21 * u13) * b22;
        u24 = (A(1,3) - b21 * u14) * b22;
        u25 = (A(1,4) - b21 * u15) * b22;
        b31 = A(2,0);
        b32 = A(2,1) - b31 * u12;
        b33 = 1 / (A(2,2) - b31 * u13 - b32 * u23);
        u34 = (A(2,3) - b31 * u14 - b32 * u24) * b33;
        u35 = (A(2,4) - b31 * u15 - b32 * u25) * b33;
        b41 = A(3,0);
        b42 = A(3,1) - b41 * u12;
        b43 = A(3,2) - b41 * u13 - b42 * u23;
        b44 = 1 / (A(3,3) - b41 * u14 - b42 * u24 - b43 * u34);
        u45 = (A(3,4) - b41 * u15 - b42 * u25 - b43 * u35) * b44;
        b51 = A(4,0);
        b52 = A(4,1) - b51 * u12;
        b53 = A(4,2) - b51 * u13 - b52 * u23;
        b54 = A(4,3) - b51 * u14 - b52 * u24 - b53 * u34;
        b55 = 1 / (A(4,4) - b51 * u15 - b52 * u25 - b53 * u35 - b54 * u45);
        //
        d1 = f(0) * b11;
        d2 = (f(1) - b21 * d1) * b22;
        d3 = (f(2) - b31 * d1 -b32 * d2) *b33;
        d4 = (f(3) - b41 * d1 -b42 * d2 - b43 * d3) * b44;
        d5 = (f(4) - b51 * d1 -b52 * d2 - b53 * d3 - b54 * d4) * b55;
        //
        Vector5 x;
        x(4) = d5;
        x(3) = d4 - u45 * d5;
        x(2) = d3 - u34 * x(3) - u35 * d5;
        x(1) = d2 - u23 * x(2) - u24 * x(3) - u25 * d5;
        x(0) = d1 - u12 * x(1) - u13 * x(2) - u14 * x(3) - u15 * d5;

        /*double sum = 0.;
        for (int i=0; i<NVAR; ++i)
        {
            sum += x[i] * A(0,i);
        }
        if (std::abs(sum - f(0)) > 1e-9)
        {
            std::cout << "diff: " << std::abs(sum - f(0)) << std::endl;
            std::cout << "f: " << f(0) << std::endl;
            std::cout << "f: " << f(1) << std::endl;
            std::cout << "f: " << f(2) << std::endl;
            std::cout << "f: " << f(3) << std::endl;
            std::cout << "f: " << f(4) << std::endl;
            std::cout << "A: " << A(0,0) << std::endl;
            std::cout << "A: " << A(0,1) << std::endl;
            std::cout << "A: " << A(0,2) << std::endl;
            std::cout << "A: " << A(0,3) << std::endl;
            std::cout << "A: " << A(0,4) << std::endl;

            std::cout << std::endl;

            std::cout << "A: " << A(1,0) << std::endl;
            std::cout << "A: " << A(1,1) << std::endl;
            std::cout << "A: " << A(1,2) << std::endl;
            std::cout << "A: " << A(1,3) << std::endl;
            std::cout << "A: " << A(1,4) << std::endl;

            std::cout << std::endl;

            std::cout << "A: " << A(2,0) << std::endl;
            std::cout << "A: " << A(2,1) << std::endl;
            std::cout << "A: " << A(2,2) << std::endl;
            std::cout << "A: " << A(2,3) << std::endl;
            std::cout << "A: " << A(2,4) << std::endl;

            std::cout << std::endl;

            std::cout << "A: " << A(3,0) << std::endl;
            std::cout << "A: " << A(3,1) << std::endl;
            std::cout << "A: " << A(3,2) << std::endl;
            std::cout << "A: " << A(3,3) << std::endl;
            std::cout << "A: " << A(3,4) << std::endl;

            std::cout << std::endl;

            std::cout << "A: " << A(4,0) << std::endl;
            std::cout << "A: " << A(4,1) << std::endl;
            std::cout << "A: " << A(4,2) << std::endl;
            std::cout << "A: " << A(4,3) << std::endl;
            std::cout << "A: " << A(4,4) << std::endl;
        }
        assert(std::abs(sum - f(0)) <= 1e-3);*/

        return x;
    }
}
