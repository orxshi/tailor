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
    //aero_para.u_ref = std::abs(component.u);
    aero_para.u_ref = std::abs(fs.u_);
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

	Tailor::Freestream fs;
	fs.read();

	auto solver = tailor.solver();

	double reduced_freq = component.reduced_freq;
	double u = component.u;
	double v = component.v;
	double w = component.w;
	//double u_ref = std::abs(u);
	double u_ref = fs.u_;
	double chord = component.chord;
	double aoa_o_deg = component.aoa_o_deg;
	double aoa_o = Tailor::deg_to_rad(aoa_o_deg); 

	double omega = 2. * reduced_freq * u_ref / chord; // angular frequency.

	double d_aoa_n = aoa_o * std::sin(omega * solver->nsolve() * solver->dt()); // change in pitch from mean.

	assert(solver->nsolve() != 0);

	double d_aoa_nm1 = aoa_o * std::sin(omega * (solver->nsolve() - 1) * solver->dt()); // change in pitch from mean.

	double d_aoa = d_aoa_n - d_aoa_nm1;

	//std::cout << "d_aoa: " << d_aoa << std::endl;
	//std::cout << "omega: " << omega << std::endl;
	//std::cout << "aoa_o: " << aoa_o << std::endl;
	//std::cout << "nassemble: " << tailor.assembler()->nassemble() << std::endl;
	//std::cout << "arg: " << omega * tailor.assembler()->nassemble() << std::endl;
	//std::cout << "sin: " << std::sin(omega * tailor.assembler()->nassemble() * 1) << std::endl;
	//std::cout << "d_aoa: " << -d_aoa << " " << Tailor::rad_to_deg(-d_aoa) << " " << -Tailor::rad_to_deg(aoa_mean + d_aoa_n) << std::endl;
	//std::cout << "aoa_total_deg: " << aoa_mean_deg + Tailor::rad_to_deg(d_aoa) << std::endl;

	int axis = component.rotaxis_;

	tailor.rotate(Tailor::Tag(1), -d_aoa, axis, component.pivot_);

	//if (tailor.comm().rank() == 0)
	//{
	//std::ofstream of;
	//of.open("angle.dat", std::ios_base::app);

	//of << Tailor::rad_to_deg(d_aoa) << " " << Tailor::rad_to_deg(aoa_mean + d_aoa_n) << std::endl;
	//of.close();
	//}
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&rotate, &set_aero_para);

    return 0;
}
