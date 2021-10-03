#include "tailor.h"

std::vector<Tailor::AeroCoefPara> set_aero_para()
{
    Tailor::Freestream fs;
    fs.read();

    Tailor::Component component;
    component.read(Tailor::Tag(1));

    Tailor::AeroCoefPara aero_para;
    aero_para.p_ref = fs.p_;
    aero_para.rho_ref = fs.rho_;
    aero_para.u_ref = fs.u_;
    aero_para.area_ref = 1.;
    aero_para.moment_center = component.pivot_;
    aero_para.moment_length = component.chord;

    std::vector<Tailor::AeroCoefPara> aero_para_vec(2, aero_para);

    return aero_para_vec;
}

void rotate(Tailor::Tailor& tailor)
{
    Tailor::Component component;
    component.read(Tailor::Tag(1));

    auto solver = tailor.solver();

    double reduced_freq = component.reduced_freq;
    double u = component.u;
    double v = component.v;
    double w = component.w;
    Tailor::Vector3 U(u, v, w);
    double u_ref = U.len();
    double chord = component.chord;
    double aoa_mean_deg = component.aoa_mean_deg; 
    double aoa_o_deg = component.aoa_o_deg;
    double aoa_mean = Tailor::deg_to_rad(aoa_mean_deg); 
    double aoa_o = Tailor::deg_to_rad(aoa_o_deg); 

    double om = 2. * reduced_freq * u_ref / chord; // rad/s

    double aoa = aoa_mean + aoa_o * std::sin(om * solver->nsolve() * solver->dt());

    int axis = component.rotaxis_;

    tailor.rotate(Tailor::Tag(1), aoa, axis, component.pivot_);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&rotate, &set_aero_para);

    return 0;
}
