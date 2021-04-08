#ifndef CORE_H
#define	CORE_H

//#define TAILOR_N_DIM 3
//#define TAILOR_BIG_POS_NUM 1e9
//#define TAILOR_BIG_NEG_NUM -1e9
//#define MASTER 0
//#define TAILOR_ZERO 1e-12
//#define ARMA_DONT_USE_WRAPPER

namespace Tailor
{
    const double TAILOR_ZERO = 1e-15;
    const int TAILOR_N_DIM = 3;
    const double TAILOR_BIG_POS_NUM = 1e9;
    const double TAILOR_BIG_NEG_NUM = -1e9;
    //static int MASTER = 0;
    const int NVAR = 5;
    //const double GAMMA = 1.4; // ratio of specific heats
    const double PI = 3.14159265359;

    enum class motion_t
    {
        translation = 0,
        rotation = 1
    };

    enum class HoleCutType
    {
        undefined,
        holemap,
        direct,
        implicit,
    };

    //enum class LoadEstimType
    //{
        //area = 0,
        //hybrid = 1,
        //solver = 2,
        //minmesh = 3,
    //};

    enum OGA_cell_type_t
    {
        undefined = -1,
        non_resident = 0,
        receptor = 1,
        field = 2,
        hole = 3,
        mandat_receptor = 4,
        orphan = 5,
        hole_candidate = 6,
        ghost = 7,
    };

    /*enum class boundary_t
    {
        undefined = -1,
        wall = 1,
        dirichlet = 2,
        empty = 3,
        interior = 4,
        partition = 4,
    };*/

    static double deg_to_rad(double deg)
    {
        return deg * PI / 180.;
    }

    static double rad_to_deg(double rad)
    {
        return rad / PI * 180.;
    }
}

#endif
