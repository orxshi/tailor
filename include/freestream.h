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
        //double aoa_air_x_;
        //double aoa_air_z_;
        //double rhoinf_;
        //double pinf_;
        //double machair_;
        //double velair_;
        double gamma_;
        double rho_;
        double p_;
        double u_;
        double v_;
        double w_;

        void read();

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            //ar & aoa_air_x_;
            //ar & aoa_air_z_;
            ar & rho_;
            ar & p_;
            //ar & machair_;
            //ar & velair_;
            ar & gamma_;
            ar & u_;
            ar & v_;
            ar & w_;
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
        double x;
        double ul;
        double ur;
        double pl;
        double pr;
        double rhol;
        double rhor;

        void read(const Tag& meshtag);
    };

    struct Component
    {
        bool rotation_;
        bool oscillation;
        int rotaxis_;
        int rpm_;
        double mach_;
        double dirx_;
        double dirz_;
        Vector3 pivot_;
        double u;
        double v;
        double w;
        double reduced_freq;
        double chord;
        double aoa_mean_deg;
        double aoa_o_deg;
        bool compute_pres_coef;
        bool compute_force_coef;
        bool compute_moment_coef;
        int priority;
        double aoa;

        void read(const Tag& mtag);
    };

    char* cstr(std::string sdesc, std::string sub);
}

#endif
