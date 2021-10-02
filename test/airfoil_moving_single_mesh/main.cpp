#include "tailor.h"

std::vector<Tailor::AeroCoefPara> set_aero_para()
{
    Tailor::Freestream fs;
    fs.read();

    Tailor::AeroCoefPara aero_para;
    aero_para.p_ref = fs.p_;
    aero_para.rho_ref = fs.rho_;
    aero_para.u_ref = fs.u_;

    std::vector<Tailor::AeroCoefPara> aero_para_vec(1, aero_para);

    return aero_para_vec;
}

void rotate(Tailor::Tailor& tailor)
{
    Tailor::Component compo;
    compo.read(Tailor::Tag(0));

    double reduced_freq = compo.reduced_freq;
    double u = compo.u;
    double v = compo.v;
    double w = compo.w;
    Tailor::Vector3 U(u, v, w);
    double u_ref = U.len();
    double chord = compo.chord;
    double aoa_mean_deg = compo.aoa_mean_deg; 
    double aoa_o_deg = compo.aoa_o_deg;
    double aoa_mean = Tailor::deg_to_rad(aoa_mean_deg); 
    double aoa_o = Tailor::deg_to_rad(aoa_o_deg); 

    double om = 2. * reduced_freq * u_ref / chord; // rad/s

    double aoa = aoa_mean + aoa_o * std::sin(om * solver()->nsolve() * solver()->dt());

    int axis = compo.rotaxis_;
    Tailor::Vector3 pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));

    tailor.rotate(Tailor::Tag(1), azimuth, axis, pivot);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&rotate, &set_aero_para);

    return 0;
}
