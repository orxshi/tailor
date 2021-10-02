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

void dummy(Tailor::Tailor& tailor)
{
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&dummy, &set_aero_para);

    return 0;
}
