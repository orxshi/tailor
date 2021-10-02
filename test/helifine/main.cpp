#include "tailor.h" 

std::vector<Tailor::AeroCoefPara> set_aero_para()
{
    Tailor::Freestream fs;
    fs.read();

    double INCH_TO_M = 0.0254;

    double uncut_wing_length = 33.88 * INCH_TO_M;
    double root_cut_out = 0.24 * uncut_wing_length;
    double R = uncut_wing_length;

    double rpm = 2000.;
    double om = rpm * 2. * Tailor::PI / 60.; // rad/s

    double stip = R * om; // wing tip speed
    double A = Tailor::PI * R * R;

    double fuslen = 78.57 * INCH_TO_M / 1.997;
    Tailor::Vector3 hub(0.696 * fuslen, 0., 0.322 * fuslen);

    Tailor::AeroCoefPara aero_para;
    aero_para.p_ref = fs.p_;
    aero_para.rho_ref = fs.rho_;
    aero_para.u_ref = stip;
    aero_para.area_ref = A;
    aero_para.moment_length = R;
    aero_para.moment_center = hub;
    
    // I am using the aerodynamic parameters of the fuselage for all the components
    // since I don't need moment coefficients on the blades and the hub.
    std::vector<Tailor::AeroCoefPara> aero_para_vec(6, aero_para);

    return aero_para_vec;
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
    tailor.make(&rotate, &set_aero_para);

    return 0;
}
