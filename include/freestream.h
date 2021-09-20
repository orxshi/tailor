#ifndef TAILOR_FREESTREAM_H
#define TAILOR_FREESTREAM_H

#include <fstream>
#include "boost/program_options.hpp"
#include "matrix.h"
#include "tag.h"

namespace Tailor
{
    struct Freestream
    {
        double aoa_air_x_;
        double aoa_air_z_;
        double rhoinf_;
        double pinf_;
        double machair_;
        double velair_;
        double gamma_;

        void read();

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & aoa_air_x_;
            ar & aoa_air_z_;
            ar & rhoinf_;
            ar & pinf_;
            ar & machair_;
            ar & velair_;
            ar & gamma_;
        }
    };

    struct FlowInit
    {
        std::string type;
        double cnt_x;
        double cnt_y;
        double cnt_z;
        double strength;
        double rho;
        double p;
        double u;
        double v;
        double w;

        void read(const Tag& meshtag);
    };

    struct Component
    {
        bool rotation_;
        int rotaxis_;
        int rpm_;
        double mach_;
        double dirx_;
        double dirz_;
        Vector3 pivot_;
        double u;
        double v;
        double w;

        void read(const Tag& mtag);
    };

    char* cstr(std::string sdesc, std::string sub);
}

#endif
