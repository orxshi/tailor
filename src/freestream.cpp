#include "freestream.h"

namespace Tailor
{
    char* cstr(std::string sdesc, std::string sub)
    {
        std::string full = sdesc;
        full.append(".");
        full.append(sub);
        char* str = new char [full.length() + 1];
        str = std::strcpy(str, full.c_str());
        return str;
    }

    void Freestream::read()
    {
        namespace po = boost::program_options;
        po::options_description op;

        po::options_description desc{"freestream"};
        desc.add_options()
            //("freestream.mach-air", po::value<double>()->default_value(0), "Mach air")
            //("freestream.vel-air", po::value<double>()->default_value(0), "Vel air")
            //("freestream.aoa-air-x", po::value<double>(), "Angle of attack")
            //("freestream.aoa-air-z", po::value<double>(), "Angle of attack")
            //("freestream.rho", po::value<double>(), "")
            //("freestream.u", po::value<double>(), "")
            //("freestream.v", po::value<double>(), "")
            //("freestream.w", po::value<double>(), "")
            //("freestream.p", po::value<double>(), "")
            ("freestream.gamma", po::value<double>(), "")
            ;

        op.add(desc);
        std::ifstream settings_file("freestream.ini");

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        //rho_ = vm["freestream.rho"].as<double>();
        //p_ = vm["freestream.p"].as<double>();
        //u_ = vm["freestream.u"].as<double>();
        //v_ = vm["freestream.v"].as<double>();
        //w_ = vm["freestream.w"].as<double>();
        //machair_ = vm["freestream.mach-air"].as<double>();
        //velair_ = vm["freestream.vel-air"].as<double>();
        //aoa_air_x_ = vm["freestream.aoa-air-x"].as<double>();
        //aoa_air_z_ = vm["freestream.aoa-air-z"].as<double>();
        gamma_ = vm["freestream.gamma"].as<double>();
    }

    void Component::read(const Tag& mtag)
    {
        namespace po = boost::program_options;
        po::options_description op;

        std::string sdesc = "component ";
        sdesc.append(std::to_string(mtag()));

        po::options_description desc{sdesc};
        desc.add_options()
            (cstr(sdesc, "rotation"), po::value<bool>()->default_value(false), "")
            (cstr(sdesc, "oscillation"), po::value<bool>()->default_value(false), "")
            (cstr(sdesc, "rotaxis"), po::value<int>()->default_value(0), "")
            (cstr(sdesc, "rpm"), po::value<int>()->default_value(0), "")
            (cstr(sdesc, "pivotx"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "pivoty"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "pivotz"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "mach"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "dirx"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "dirz"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "u"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "v"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "w"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "reduced-freq"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "chord"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "aoa-mean-deg"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "aoa-o-deg"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "compute-pres-coef"), po::value<bool>()->default_value(false), "")
            (cstr(sdesc, "compute-force-coef"), po::value<bool>()->default_value(false), "")
            (cstr(sdesc, "compute-moment-coef"), po::value<bool>()->default_value(false), "")
            ;

        op.add(desc);
        std::string s = "component.ini";
        std::ifstream settings_file(s);

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        rotation_ = vm[cstr(sdesc, "rotation")].as<bool>();
        oscillation = vm[cstr(sdesc, "oscillation")].as<bool>();
        rotaxis_ = vm[cstr(sdesc, "rotaxis")].as<int>();
        rpm_ = vm[cstr(sdesc, "rpm")].as<int>();
        double pivotx = vm[cstr(sdesc, "pivotx")].as<double>();
        double pivoty = vm[cstr(sdesc, "pivoty")].as<double>();
        double pivotz = vm[cstr(sdesc, "pivotz")].as<double>();
        mach_ = vm[cstr(sdesc, "mach")].as<double>();
        dirx_ = vm[cstr(sdesc, "dirx")].as<double>();
        dirz_ = vm[cstr(sdesc, "dirz")].as<double>();
        pivot_ = Vector3(pivotx, pivoty, pivotz);
        u = vm[cstr(sdesc, "u")].as<double>();
        v = vm[cstr(sdesc, "v")].as<double>();
        w = vm[cstr(sdesc, "w")].as<double>();
        reduced_freq = vm[cstr(sdesc, "reduced-freq")].as<double>();
        chord = vm[cstr(sdesc, "chord")].as<double>();
        aoa_mean_deg = vm[cstr(sdesc, "aoa-mean-deg")].as<double>();
        aoa_o_deg = vm[cstr(sdesc, "aoa-o-deg")].as<double>();
        compute_pres_coef = vm[cstr(sdesc, "compute-pres-coef")].as<bool>();
        compute_force_coef = vm[cstr(sdesc, "compute-force-coef")].as<bool>();
        compute_moment_coef = vm[cstr(sdesc, "compute-moment-coef")].as<bool>();
    }

    void FlowInit::read(const Tag& meshtag)
    {
        namespace po = boost::program_options;
        po::options_description op;

        std::string sdesc = "mesh ";
        sdesc.append(std::to_string(meshtag()));

        po::options_description desc{sdesc};
        desc.add_options()
            (cstr(sdesc, "type"), po::value<std::string>()->required(), "")
            (cstr(sdesc, "cnt_x"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "cnt_y"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "cnt_z"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "strength"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "rho"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "p"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "u"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "v"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "w"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "ul"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "ur"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "pl"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "pr"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "rhol"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "rhor"), po::value<double>()->default_value(0), "")
            (cstr(sdesc, "x"), po::value<double>()->default_value(0), "")
            ;

        op.add(desc);
        std::string fn = "flow_init.ini";
        std::ifstream settings_file(fn);

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        type = vm[cstr(sdesc, "type")].as<std::string>();
        cnt_x = vm[cstr(sdesc, "cnt_x")].as<double>();
        cnt_y = vm[cstr(sdesc, "cnt_y")].as<double>();
        cnt_z = vm[cstr(sdesc, "cnt_z")].as<double>();
        strength = vm[cstr(sdesc, "strength")].as<double>();
        rho = vm[cstr(sdesc, "rho")].as<double>();
        p = vm[cstr(sdesc, "p")].as<double>();
        u = vm[cstr(sdesc, "u")].as<double>();
        v = vm[cstr(sdesc, "v")].as<double>();
        w = vm[cstr(sdesc, "w")].as<double>();
        ul = vm[cstr(sdesc, "ul")].as<double>();
        ur = vm[cstr(sdesc, "ur")].as<double>();
        pl = vm[cstr(sdesc, "pl")].as<double>();
        pr = vm[cstr(sdesc, "pr")].as<double>();
        rhol = vm[cstr(sdesc, "rhol")].as<double>();
        rhor = vm[cstr(sdesc, "rhor")].as<double>();
        x = vm[cstr(sdesc, "x")].as<double>();
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
