#include "freestream.h"

namespace Tailor
{
    void Freestream::read()
    {
        namespace po = boost::program_options;
        po::options_description op;

        po::options_description desc{"freestream"};
        desc.add_options()
            ("freestream.mach-air", po::value<double>()->default_value(0), "Mach air")
            ("freestream.vel-air", po::value<double>()->default_value(0), "Vel air")
            ("freestream.aoa-air-x", po::value<double>(), "Angle of attack")
            ("freestream.aoa-air-z", po::value<double>(), "Angle of attack")
            ("freestream.p", po::value<double>(), "Reference pressure")
            ("freestream.rho", po::value<double>(), "Reference density")
            ("freestream.gamma", po::value<double>(), "Reference density")
            ;

        op.add(desc);
        std::ifstream settings_file("freestream.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        rhoinf_ = vm["freestream.rho"].as<double>();
        pinf_ = vm["freestream.p"].as<double>();
        machair_ = vm["freestream.mach-air"].as<double>();
        velair_ = vm["freestream.vel-air"].as<double>();
        aoa_air_x_ = vm["freestream.aoa-air-x"].as<double>();
        aoa_air_z_ = vm["freestream.aoa-air-z"].as<double>();
        gamma_ = vm["freestream.gamma"].as<double>();
    }

    void Component::read(const Tag& mtag)
    {
        namespace po = boost::program_options;
        po::options_description op;

        po::options_description desc{"compo"};
        desc.add_options()
            ("rotation", po::value<bool>()->default_value(false), "")
            ("rotaxis", po::value<int>()->default_value(0), "")
            ("rpm", po::value<int>()->default_value(0), "")
            ("pivotx", po::value<double>()->default_value(0), "")
            ("pivoty", po::value<double>()->default_value(0), "")
            ("pivotz", po::value<double>()->default_value(0), "")
            ("mach", po::value<double>()->default_value(0), "")
            ("dirx", po::value<double>()->default_value(0), "")
            ("dirz", po::value<double>()->default_value(0), "")
            ;

        op.add(desc);
        std::string s = "compo_";
        s.append(std::to_string(mtag()));
        s.append(".ini");
        std::ifstream settings_file(s);

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        rotation_ = vm["rotation"].as<bool>();
        rotaxis_ = vm["rotaxis"].as<int>();
        rpm_ = vm["rpm"].as<int>();
        double pivotx = vm["pivotx"].as<double>();
        double pivoty = vm["pivoty"].as<double>();
        double pivotz = vm["pivotz"].as<double>();
        mach_ = vm["mach"].as<double>();
        dirx_ = vm["dirx"].as<double>();
        dirz_ = vm["dirz"].as<double>();
        pivot_ = vec3<double>(pivotx, pivoty, pivotz);
    }
}
