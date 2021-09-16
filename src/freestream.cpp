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
            ("u", po::value<double>()->default_value(0), "")
            ("v", po::value<double>()->default_value(0), "")
            ("w", po::value<double>()->default_value(0), "")
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
        pivot_ = Vector3(pivotx, pivoty, pivotz);
        u = vm["u"].as<double>();
        v = vm["v"].as<double>();
        w = vm["w"].as<double>();
    }

    void FlowInit::read(const Tag& meshtag)
    {
        auto copt = [](std::string sdesc, std::string sub)
        {
            std::string full = sdesc;
            full.append(".");
            full.append(sub);
            return full.c_str();
        };

        auto sopt = [](std::string sdesc, std::string sub)
        {
            std::string full = sdesc;
            full.append(".");
            full.append(sub);
            return full;
        };

        namespace po = boost::program_options;
        po::options_description op;

        std::string sdesc = "mesh ";
        sdesc.append(std::to_string(meshtag()));
        //auto cdesc = sdesc.c_str();
        
        po::options_description desc{sdesc};
        desc.add_options()
            (copt(sdesc, "type"), po::value<std::string>()->required(), "")
            (copt(sdesc, "cnt_x"), po::value<double>(), "")
            //(copt(sdesc, "cnt_y"), po::value<double>(), "")
            //(copt(sdesc, "cnt_z"), po::value<double>(), "")
            //(copt(sdesc, "strength"), po::value<double>(), "")
            //(copt(sdesc, "rho"), po::value<double>(), "")
            //(copt(sdesc, "p"), po::value<double>(), "")
            //(copt(sdesc, "u"), po::value<double>(), "")
            //(copt(sdesc, "v"), po::value<double>(), "")
            //(copt(sdesc, "w"), po::value<double>(), "")
            ;

        op.add(desc);
        std::string fn = "flow_init.ini";
        //fn.append(std::to_string(meshtag()));
        //fn.append(".ini");
        std::ifstream settings_file(fn);

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        type = vm[sopt(sdesc, "type")].as<std::string>();
        cnt_x = vm[sopt(sdesc, "cnt_x")].as<double>();
        //cnt_y = vm[sopt(sdesc, "cnt_y")].as<double>();
        //cnt_z = vm[sopt(sdesc, "cnt_z")].as<double>();
        //strength = vm[sopt(sdesc, "strength")].as<double>();
        //rho = vm[sopt(sdesc, "rho")].as<double>();
        //p = vm[sopt(sdesc, "p")].as<double>();
        //u = vm[sopt(sdesc, "u")].as<double>();
        //v = vm[sopt(sdesc, "v")].as<double>();
        //w = vm[sopt(sdesc, "w")].as<double>();
    }

    //void GaussianInit::read()
    //{
    //    namespace po = boost::program_options;
    //    po::options_description op;

    //    po::options_description desc{"ginit"};
    //    desc.add_options()
    //        ("ginit.cnt_x", po::value<double>()->required(), "")
    //        ("ginit.cnt_y", po::value<double>()->required(), "")
    //        ("ginit.cnt_z", po::value<double>()->required(), "")
    //        ("ginit.strength", po::value<double>()->required(), "")
    //        //("ginit.p", po::value<double>(), "")
    //        //("ginit.u", po::value<double>(), "")
    //        //("ginit.v", po::value<double>(), "")
    //        //("ginit.w", po::value<double>(), "")
    //        //("ginit.rho", po::value<double>(), "")
    //        ;

    //    op.add(desc);
    //    std::ifstream settings_file("gaussian_init.ini");

    //    boost::program_options::variables_map vm;
    //    po::store(po::parse_config_file(settings_file, op, true), vm);
    //    po::notify(vm);

    //    cnt_x = vm["ginit.cnt_x"].as<double>();
    //    cnt_y = vm["ginit.cnt_y"].as<double>();
    //    cnt_z = vm["ginit.cnt_z"].as<double>();
    //    strength = vm["ginit.strength"].as<double>();
    //    //p = vm["freestream.rho"].as<double>();
    //    //u = vm["freestream.rho"].as<double>();
    //    //v = vm["freestream.rho"].as<double>();
    //    //w = vm["freestream.rho"].as<double>();
    //    //rho = vm["freestream.p"].as<double>();
    //}
}
