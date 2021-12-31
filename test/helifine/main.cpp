#include "tailor.h" 

using namespace Tailor;

double current_azimuth_1 = 0.5 * PI;
double current_azimuth_2 = 1.0 * PI;
double current_azimuth_3 = 1.5 * PI;
double current_azimuth_4 = 0.0 * PI;

double old_cyclic_pitch_1 = 0.;
double old_cyclic_pitch_2 = 0.;
double old_cyclic_pitch_3 = 0.;
double old_cyclic_pitch_4 = 0.;

Vector3 blade_vector_1 = Vector3(0., 1., 0.); 
Vector3 blade_vector_2 = Vector3(-1., 0., 0.); 
Vector3 blade_vector_3 = Vector3(0., -1., 0.); 
Vector3 blade_vector_4 = Vector3(1., 0., 0.); 

std::vector<AeroCoefPara> set_aero_para()
{
    Freestream fs;
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
    Vector3 hub(0.696 * fuslen, 0., 0.322 * fuslen);

    AeroCoefPara aero_para;
    aero_para.p_ref = fs.p_;
    aero_para.rho_ref = fs.rho_;
    aero_para.u_ref = stip;
    aero_para.area_ref = A;
    aero_para.moment_length = R;
    aero_para.moment_center = hub;
    
    // I am using the aerodynamic parameters of the fuselage for all the components
    // since I don't need moment coefficients on the blades and the hub.
    std::vector<AeroCoefPara> aero_para_vec(6, aero_para);

    return aero_para_vec;
}

void rotate(Tailor::Tailor& tailor)
{
    RotationMatrix rm;

    Component compo;
    compo.read(Tag(1));

    double rpm = compo.rpm_;
    double om = rpm * 2. * PI / 60.; // rad/s
    double delta_azimuth = om * tailor.solver()->dt(); // rad/s * time step
    int axis = compo.rotaxis_;
    Vector3 pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));
    Vector3 origin(0., 0., 0.);
    Vector3 z(0., 0., 1.);

    // rotate blade vectors azimuthally.
    current_azimuth_1 += delta_azimuth;
    current_azimuth_2 += delta_azimuth;
    current_azimuth_3 += delta_azimuth;
    current_azimuth_4 += delta_azimuth;

    rm.rotate(delta_azimuth, 2, pivot, blade_vector_1);
    rm.rotate(delta_azimuth, 2, pivot, blade_vector_2);
    rm.rotate(delta_azimuth, 2, pivot, blade_vector_3);
    rm.rotate(delta_azimuth, 2, pivot, blade_vector_4);

    // calculate cyclic pitch
    double cyclic_pitch_1 = -0.1 * std::sin(current_azimuth_1) + 0.2 * std::cos(current_azimuth_1);
    double cyclic_pitch_2 = -0.1 * std::sin(current_azimuth_2) + 0.2 * std::cos(current_azimuth_2);
    double cyclic_pitch_3 = -0.1 * std::sin(current_azimuth_3) + 0.2 * std::cos(current_azimuth_3);
    double cyclic_pitch_4 = -0.1 * std::sin(current_azimuth_4) + 0.2 * std::cos(current_azimuth_4);

    // pitch blade meshes cyclicly.
    tailor.rotate(Tag(1), cyclic_pitch_1 - old_cyclic_pitch_1, blade_vector_1, origin);
    tailor.rotate(Tag(2), cyclic_pitch_2 - old_cyclic_pitch_2, blade_vector_2, origin);
    tailor.rotate(Tag(3), cyclic_pitch_3 - old_cyclic_pitch_3, blade_vector_3, origin);
    tailor.rotate(Tag(4), cyclic_pitch_4 - old_cyclic_pitch_4, blade_vector_4, origin);

    // yaw blade meshes. 
    tailor.rotate(Tag(1), delta_azimuth, z, pivot);
    tailor.rotate(Tag(2), delta_azimuth, z, pivot);
    tailor.rotate(Tag(3), delta_azimuth, z, pivot);
    tailor.rotate(Tag(4), delta_azimuth, z, pivot);

    old_cyclic_pitch_1 = cyclic_pitch_1;
    old_cyclic_pitch_2 = cyclic_pitch_2;
    old_cyclic_pitch_3 = cyclic_pitch_3;
    old_cyclic_pitch_4 = cyclic_pitch_4;
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&rotate, &set_aero_para);

    return 0;
}
