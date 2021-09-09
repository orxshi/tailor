#include "tailor.h" 

Tailor::AeroCoefPara set_aero_para()
{
    Tailor::Freestream fs;
    fs.read();

    double cinf = std::sqrt(fs.gamma_ * fs.pinf_ / fs.rhoinf_);

    auto vinf_air = Tailor::Vector3(
            fs.machair_ * cinf * std::cos(Tailor::deg_to_rad(fs.aoa_air_x_)),
            fs.machair_ * cinf * std::cos(Tailor::deg_to_rad(90. - fs.aoa_air_x_)),
            fs.machair_ * cinf * std::cos(Tailor::deg_to_rad(fs.aoa_air_z_)));

    Tailor::AeroCoefPara aero_para;
    aero_para.p_ref = fs.pinf_;
    aero_para.rho_ref = fs.rhoinf_;
    aero_para.u_ref = vinf_air.len();

    return aero_para;
}

void dummy(Tailor::Tailor& tailor)
{
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(dummy);

    auto aero_para_airfoil = set_aero_para();
    std::vector<Tailor::AeroCoefPara> aero_para(1, aero_para_airfoil);
    tailor.get_aero_coef(aero_para);

    return 0;
}
