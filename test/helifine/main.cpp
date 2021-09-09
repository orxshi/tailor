#include "tailor.h" 

AeroCoefPara set_aero_para()
{
    Tailor::Freestream fs;
    fs.read();

    double INCH_TO_M = 0.0254;

    double uncut_wing_length = 33.88 * INCH_TO_M;
    double root_cut_out = 0.24 * uncut_wing_length;
    double R = uncut_wing_length;

    double rpm = 2000.;
    double om = rpm * 2. * PI / 60.; // rad/s

    double stip = R * om; // wing tip speed
    double A = PI * R * R;

    double fuslen = 78.57 * INCH_TO_M / 1.997;
    Tailor::Vector3 hub(0.696 * fuslen, 0., 0.322 * fuslen);

    AeroCoefPara aero_para;
    aero_para.p_ref = fs.pinf_;
    aero_para.rho_ref = fs.rhoinf_;
    aero_para.u_ref = stip;
    aero_para.area_ref = A;
    aero_para.moment_length = R;
    aero_para.moment_center = hub;

    return aero_para;
}

void rotate(Tailor::Tailor& tailor)
{
    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    double rpm = compo.rpm_;
    double om = rpm * 2. * Tailor::PI / 60.; // rad/s
    double azimuth = om * tailor.solver()->dt(); // rad/s * time step
    int axis = compo.rotaxis_;
    Tailor::Vector3 pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));

    tailor.rotate(Tailor::Tag(1), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(2), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(3), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(4), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(5), azimuth, axis, pivot);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(rotate);

    // I am using the aerodynamic parameters of the fuselage for all the components
    // since I don't need moment coefficients on the blades and the hub.
    auto aero_para_fuselage = set_aero_para();
    std::vector<Tailor::AeroCoefPara> aero_para(aero_para_fuselage, 6);
    tailor.get_aero_coef(aero_para);

    return 0;
}
