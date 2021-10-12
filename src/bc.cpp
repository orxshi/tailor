#include "bc.h"

namespace Tailor
{
    Matrix5 make_rot_matrix(const Vector3& n, bool verbose)
    {
        //auto l = tangent_vector(n);
        //auto m = cross(n, l);

        //double nx = n(0);
        //double ny = n(1);
        //double nz = n(2);
        //
        //NormTang nt(n);

        //double nx = n(0);
        //double ny = n(1);
        //double nz = n(2);

        //double lx = nt.l(0);
        //double ly = nt.l(1);
        //double lz = nt.l(2);

        //double mx = nt.m(0);
        //double my = nt.m(1);
        //double mz = nt.m(2);

        //double thetaz = atan2(ny, nx); // limit thetaz to first quad. otherwise consider case where n = (-1, 0, 0) which will be rotated back to (1,0,0) if we dont't limit thetaz.
        //double thetay = atan2(std::abs(nz), std::abs(nx));

        //double cosy = std::cos(thetay);
        //double cosz = std::cos(thetaz);
        //double siny = std::sin(thetay);
        //double sinz = std::sin(thetaz);

        //double lx = -sinz * cosy;
        //double ly = cosy * sinz;
        //double lz = -siny * sinz;

        //double mx = -siny;
        //double my = 0.;
        //double mz = cosy;

        //Vector3 l(lx, ly, lz);
        //Vector3 m(mx, my, mz);

        //if (verbose)
        //{
            //std::cout << "cosy * cosz: " << cosy * cosz << std::endl;
            //std::cout << "cosy * sinz: " << cosy * sinz << std::endl;
            //std::cout << "siny * cosz: " << siny * cosz << std::endl;

            //std::cout << "thetay: " << thetay << std::endl;
            //std::cout << "thetaz: " << thetaz << std::endl;

            //std::cout << "nx: " << n(0) << std::endl;
            //std::cout << "ny: " << n(1) << std::endl;
            //std::cout << "nz: " << n(2) << std::endl;

            //std::cout << "lx: " << l(0) << std::endl;
            //std::cout << "ly: " << l(1) << std::endl;
            //std::cout << "lz: " << l(2) << std::endl;

            //std::cout << "mx: " << m(0) << std::endl;
            //std::cout << "my: " << m(1) << std::endl;
            //std::cout << "mz: " << m(2) << std::endl;

            //std::cout << "dp: " << dotp(n, l) << std::endl;
            //std::cout << "dp: " << dotp(n, m) << std::endl;

            //std::cout << "len n: " << n.len() << std::endl;
            //std::cout << "len l: " << l.len() << std::endl;
            //std::cout << "len m: " << m.len() << std::endl;

            //std::cout << "cross(0): " << cross(n, l)(0) << std::endl;
            //std::cout << "cross(1): " << cross(n, l)(1) << std::endl;
            //std::cout << "cross(2): " << cross(n, l)(2) << std::endl;

            //std::cout << "cross(0): " << cross(m, n)(0) << std::endl;
            //std::cout << "cross(1): " << cross(m, n)(1) << std::endl;
            //std::cout << "cross(2): " << cross(m, n)(2) << std::endl;
        //}

        /*Matrix5 T;

        T(0,0) = 1.;
        T(0,1) = 0.;
        T(0,2) = 0.;
        T(0,3) = 0.;
        T(0,4) = 0.;

        T(1,0) = 0.;
        T(1,1) = n(0);
        T(1,2) = n(1);
        T(1,3) = n(2);
        T(1,4) = 0.;

        T(2,0) = 0.;
        T(2,1) = l(0);
        T(2,2) = l(1);
        T(2,3) = l(2);
        T(2,4) = 0.;

        T(3,0) = 0.;
        T(3,1) = m(0);
        T(3,2) = m(1); 
        T(3,3) = m(2);
        T(3,4) = 0.;

        T(4,0) = 0.;
        T(4,1) = 0.;
        T(4,2) = 0.;
        T(4,3) = 0.;
        T(4,4) = 1.;

        return T;*/


        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        double thetaz = atan2(ny, nx); // limit thetaz to first quad. otherwise consider case where n = (-1, 0, 0) which will be rotated back to (1,0,0) if we dont't limit thetaz.

        double cosz = std::cos(thetaz);
        double sinz = std::sin(thetaz);
        
        //double thetay = atan2(nz, std::abs(nx));
        double thetay = asin(nz);
        //double thetay = atan2(nz * cosz, nx);
        //double thetay = atan2(nz * cosz, std::abs(nx));

        double cosy = std::cos(thetay);
        double siny = std::sin(thetay);

        Matrix5 T;

        T(0,0) = 1.;
        T(0,1) = 0.;
        T(0,2) = 0.;
        T(0,3) = 0.;
        T(0,4) = 0.;

        T(1,0) = 0.;
        T(1,1) = cosy * cosz;
        T(1,2) = cosy * sinz;
        T(1,3) = siny;
        T(1,4) = 0.;

        T(2,0) = 0.;
        T(2,1) = -sinz;
        T(2,2) = cosz;
        T(2,3) = 0.;
        T(2,4) = 0.;

        T(3,0) = 0.;
        T(3,1) = -siny * cosz;
        T(3,2) = -siny * sinz;
        T(3,3) = cosy;
        T(3,4) = 0.;

        T(4,0) = 0.;
        T(4,1) = 0.;
        T(4,2) = 0.;
        T(4,3) = 0.;
        T(4,4) = 1.;

        //std::cout << nx << " " << ny << " " << nz <<  " " << thetaz << " " << cosz << " " << cosy << " " << T(1,1) << std::endl;
        //std::cout << "nx: " << nx << std::endl;
        //std::cout << "ny: " << ny << std::endl;
        //std::cout << "nz: " << nz << std::endl;
        //std::cout << "thetaz: " << thetaz << std::endl;
        //std::cout << "thetay: " << thetay << std::endl;

        return T;
    }
    BoundaryCondition::BoundaryCondition(): first_(true)
    {
    }

    void BoundaryCondition::farfield5(Mesh& mesh)
    {
        // toro p. 224

        for (MeshCell& mc: mesh.farfield_boundaries_)
        {
            const MeshCell& nei = mesh.cell(mc.interior_boundary());

            mc.prim_ = nei.prim();
            mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);
        }
    }

    void BoundaryCondition::farfield4(Mesh& mesh)
    {
        // https://www.researchgate.net/publication/264044118_A_Guide_to_the_Implementation_of_Boundary_Conditions_in_Compact_High-Order_Methods_for_Compressible_Aerodynamics

        double pinf = fs_.p_;
        double rhoinf = fs_.rho_;
        //double aoa_air_x = fs_.aoa_air_x_;
        //double aoa_air_z = fs_.aoa_air_z_;
        //double aoa_foil_x = fs_.aoa_foil_x_;
        //double aoa_foil_z = fs_.aoa_foil_z_;
        //double machair = fs_.machair_;
        double gamma = fs_.gamma_;
        //double machfoil = fs_.machfoil_;
        double u = fs_.u_;
        double v = fs_.v_;
        double w = fs_.w_;

        double cinf = std::sqrt(gamma * pinf / rhoinf);

        Vector3 vinf_air;
        Vector3 compo;

        //compo(0) = std::cos(deg_to_rad(fs_.aoa_air_x_));
        //compo(1) = std::cos(deg_to_rad(90. - fs_.aoa_air_x_));
        //compo(2) = std::cos(deg_to_rad(fs_.aoa_air_z_));

        vinf_air(0) = u;
        vinf_air(1) = v;
        vinf_air(2) = w;

        //if (fs_.velair_ != 0.)
        {
            //vinf_air = compo * fs_.velair_;
        }
        //else
        //{
            //vinf_air = compo * fs_.machair_ * cinf;
        //}

        for (MeshCell& mc: mesh.farfield_boundaries_)
        {
            assert(mc.face().size() == 1);

            const MeshCell& nei = mesh.cell(mc.interior_boundary());
            auto mf = std::find_if(nei.face().begin(), nei.face().end(), [&](const MeshFace& f){return f.tag() == mc.face()[0].tag();});
            assert(mf != nei.face().end());
            auto n = mf->face().normal();

            //const auto& mf = mc.face()[0];
            //Vector3 n = mf.face().normal();

            //Vector3 vinf = vinf_air - normalize(n) * mesh.face(mf.tag()).vgn();
            //Vector3 vinf = vinf_air - normalize(n) * mf->vgn();
            //Vector3 vinf = vinf_air - mf->vgn();
            //Vector3 vinf = vinf_air - (n * mf->vgn());
            //Vector3 vinf = vinf_air - mf->vf();
            Vector3 vinf = vinf_air;
            //if (mf->vf()(0) != 0. || vinf(1) != 0.)
            //{
            //    std::cout << "vf(0): " << mf->vf()(0) << std::endl;
            //    std::cout << "vinf(0): " << vinf(0) << std::endl;
            //    std::cout << "vinf(1): " << vinf(1) << std::endl;
            //    std::cout << "vinf(2): " << vinf(2) << std::endl;
            //}
            //assert(mf->vf()(0) == 0.);
            //assert(mf->vf()(1) == 0.);
            //assert(mf->vf()(2) == 0.);

            //assert(std::abs(vinf(1)) < TAILOR_ZERO);
            //assert(std::abs(vinf(2)) < TAILOR_ZERO);

            const MeshCell& interior = mesh.cell(mc.interior_boundary());

            if (interior.oga_cell_type() == OGA_cell_type_t::mandat_receptor || interior.oga_cell_type() == OGA_cell_type_t::non_resident)
            {
                /*if (mesh.tag()() == 0)
                {
                    if (mc.tag()() == 29)
                    {
                        std::cout << "interior: " << interior.tag()() << std::endl;
                        std::cout << "type: " << static_cast<int>(interior.oga_cell_type()) << std::endl;
                        mesh.print_as_vtk("asu.vtk");
                        assert(false);
                    }
                }*/
                continue;
            }

            double rhoi= interior.prim(0);
            double pi = interior.prim(4);
            if (pi < 0. || std::isnan(pi))
            {
                std::cout << "prim(0): " << interior.prim(0) << std::endl;
                std::cout << "prim(1): " << interior.prim(1) << std::endl;
                std::cout << "prim(2): " << interior.prim(2) << std::endl;
                std::cout << "prim(3): " << interior.prim(3) << std::endl;
                std::cout << "prim(4): " << interior.prim(4) << std::endl;
            }
            assert(pi >= 0.);
            Vector3 vi(interior.prim(1), interior.prim(2), interior.prim(3));
            double vin = dot(vi, n);
            double vinfn = dot(vinf, n);
            double ci = std::sqrt(gamma * pi / rhoi);

            double rplus, rminus;

            bool subsonic = false;
            if ((std::abs(vin) / ci) < 1) {
                subsonic = true;
            }

            if (vin < 0.) // inflow
            {
                rminus = vinfn - 2. * cinf / (gamma - 1.);
                if (subsonic)
                {
                    rplus = vin + 2. * ci / (gamma - 1.);
                }
                else
                {
                    rplus = vinfn + 2. * cinf / (gamma - 1.);
                }

                double vb = 0.5 * (rplus + rminus);
                double cb = (gamma - 1.) * (rplus - rminus) / 4.;

                Vector3 v = vinf + n * (vb - vinfn);

                //if (std::abs(v(1)) > TAILOR_ZERO)
                //if (interior.tag()() == 23475)
                //{
                //    std::cout << "inflow" << std::endl;
                //    std::cout << "subsonic: " << subsonic << std::endl;
                //    std::cout << "n(0): " << n(0) << std::endl;
                //    std::cout << "n(1): " << n(1) << std::endl;
                //    std::cout << "n(2): " << n(2) << std::endl;
                //    std::cout << "v(1): " << v(1) << std::endl;
                //    std::cout << "vinf(1): " << vinf(1) << std::endl;
                //    std::cout << "vb: " << vb << std::endl;
                //    std::cout << "vi(0): " << vi(0) << std::endl;
                //    std::cout << "vi(1): " << vi(1) << std::endl;
                //    std::cout << "vi(2): " << vi(2) << std::endl;
                //    assert(false);
                //}
                //assert(std::abs(v(1)) < TAILOR_ZERO);
                //assert(std::abs(v(2)) < TAILOR_ZERO);

                double sb = std::pow(cinf, 2.) / gamma / std::pow(rhoinf, gamma - 1.);
                double rhob = std::pow(cb, 2.) / gamma / sb;
                double pb = rhob * std::pow(cb, 2.) / gamma; 

                if (pb < 0. || isnan(pb))
                {
                    std::cout << "rhoinf: " << rhoinf << std::endl;
                    std::cout << "gamma: " << gamma << std::endl;
                    std::cout << "cinf: " << cinf << std::endl;
                    std::cout << "sb: " << sb << std::endl;
                    std::cout << "rhob: " << sb << std::endl;
                    std::cout << "pb: " << pb << std::endl;
                }
                assert(pb >= 0.);
                assert(rhob > 0.);

                assert(rhob != 0.);
                mc.prim_(0) = rhob;
                mc.prim_(1) = v(0);
                mc.prim_(2) = v(1);
                mc.prim_(3) = v(2);
                mc.prim_(4) = pb;

                assert(mc.prim_(0) == rhoinf);
                assert(mc.prim_(1) == u);
                assert(mc.prim_(2) == 0.);
                assert(mc.prim_(3) == w);
                assert(mc.prim_(4) == pinf);

                mc.cons_sp1_ = prim_to_cons(mc.prim_, gamma);
            }
            else // outflow
            {
                rplus = vin + 2. * ci / (gamma - 1.);
                if (subsonic)
                {
                    rminus = vinfn - 2. * cinf / (gamma - 1.);
                }
                else
                {
                    rminus = vin - 2. * ci / (gamma - 1.);
                }

                double vb = 0.5 * (rplus + rminus);
                double cb = (gamma - 1.) * (rplus - rminus) / 4.;

                Vector3 v = vi + n * (vb - vin);


                //if (std::abs(v(1)) > TAILOR_ZERO)
                //if (interior.tag()() == 23475)
                //{
                //    std::cout << "outflow" << std::endl;
                //    std::cout << "subsonic: " << subsonic << std::endl;
                //    std::cout << "n(0): " << n(0) << std::endl;
                //    std::cout << "n(1): " << n(1) << std::endl;
                //    std::cout << "n(2): " << n(2) << std::endl;
                //    std::cout << "v(1): " << v(1) << std::endl;
                //    std::cout << "vinf(1): " << vinf(1) << std::endl;
                //    std::cout << "vb: " << vb << std::endl;
                //    std::cout << "vi(0): " << vi(0) << std::endl;
                //    std::cout << "vi(1): " << vi(1) << std::endl;
                //    std::cout << "vi(2): " << vi(2) << std::endl;
                //    assert(false);
                //}
                //assert(std::abs(v(1)) < TAILOR_ZERO);
                //assert(std::abs(v(2)) < TAILOR_ZERO);

                double sb = std::pow(ci, 2.) / gamma / std::pow(rhoi, gamma - 1.);
                double rhob = std::pow(cb, 2.) / gamma / sb;
                double pb = rhob * std::pow(cb, 2.) / gamma; 

                //if (pb < 0. || std::isnan(pb))
                //{
                //    std::cout << "pb: " << pb << std::endl;
                //    std::cout << "rhob: " << rhob << std::endl;
                //    std::cout << "cb: " << cb << std::endl;
                //    std::cout << "ci: " << ci << std::endl;
                //    std::cout << "vin: " << vin << std::endl;
                //    std::cout << "vinfn: " << vinfn << std::endl;
                //    std::cout << "gamma: " << gamma << std::endl;
                //    std::cout << "mf.vgn(): " << mf->vgn() << std::endl;
                //}
                assert(pb >= 0.);
                assert(rhob > 0.);

                assert(rhob != 0.);
                mc.prim_(0) = rhob;
                mc.prim_(1) = v(0);
                mc.prim_(2) = v(1);
                mc.prim_(3) = v(2);
                mc.prim_(4) = pb;

                assert(mc.prim_(0) == rhoinf);
                if (mc.prim_(1) != u)
                {
                    std::cout << "vi: " << vi(0) << " " << vi(1) << " " << vi(2) << std::endl;
                    std::cout << "vb: " << vb << std::endl;
                    std::cout << "vin: " << vin << std::endl;
                }
                assert(mc.prim_(1) == u);
                assert(mc.prim_(2) == 0.);
                assert(mc.prim_(3) == w);
                assert(mc.prim_(4) == pinf);

                mc.cons_sp1_ = prim_to_cons(mc.prim_, gamma);
            }
        }
    }

    /*void BoundaryCondition::farfield3(Mesh& mesh)
    {
        // https://github.com/su2code/SU2_EDU/blob/master/src/solver_direct_mean.cpp#L2384

        //double vb = 0.;
        std::ifstream in;
        in.open("init_farfield.bc");
        if (!in.is_open()) {
            return;
        }

        double rhoinf;
        double pinf;
        double machinf;
        double aoa;

        in >> rhoinf;
        in >> pinf;
        in >> machinf;
        in >> aoa;

        in.close();
        
        double cinf = std::sqrt(GAMMA * pinf / rhoinf);
        //double vinf = machinf * cinf;
        Vector3 vinf(machinf * cinf * std::cos(deg_to_rad(aoa)), machinf * cinf * std::sin(deg_to_rad(aoa)), 0.);
        double sinf = std::pow(rhoinf, GAMMA) / pinf;

        for (MeshCell& mc: mesh.farfield_boundaries_)
        {
            assert(mc.face().size() == 1);
            Vector3 normal = mc.face()[0].face().normal();

            //Vector3 vinf1 = normal * vinf;
            Vector3 vinf1 = normal * dotp(vinf, normal);

            const MeshCell& interior = mesh.cell(mc.interior_boundary());

            double rhobound = interior.prim(0);
            double pbound = interior.prim(4);
            assert(pbound >= 0.); // lets return u_inf values if p is negative
            Vector3 vbound1(interior.prim(1), interior.prim(2), interior.prim(3));
            double vbound = dotp(vbound1, normal);
            double cbound = std::sqrt(GAMMA * pbound / rhobound);
            double sbound = std::pow(rhobound, GAMMA) / pbound;

            //if (vn < 0.) --> inflow
                //if (vn > -c) --> subsonic
                //if (vn < -c) --> supersonic
            //if (vn > 0.) --> outflow
                //if (vn < c) --> subsonic
                //if (vn > c) --> supersonic

            double rplus, rminus, s, rho, p, v;
            //std::cout << "vinf1(0): " << vinf1(0) << std::endl;
            //std::cout << "vinf1(1): " << vinf1(1) << std::endl;
            //std::cout << "vinf1(2): " << vinf1(2) << std::endl;
            //std::cout << "vinflen: " << vinf.len() << std::endl;
            //std::cout << "normal(0): " << normal(0) << std::endl;
            //std::cout << "normal(1): " << normal(1) << std::endl;
            //std::cout << "normal(2): " << normal(2) << std::endl;
            //std::cout << "dotp(vinf, normal): " << dotp(vinf, normal) << std::endl;
            bool inflow = false;
            bool subsonic = false;

            if (dotp(vinf, normal) < 0.) // inflow
            {
                inflow = true;
                if (dotp(vinf, normal) > -cinf) // subsonic
                {
                    subsonic = true;
                    //std::cout << "subsonic inflow" << std::endl;
                    rminus = vbound - 2. * cbound / (GAMMA - 1.);

                    double e = spec_inter_energy(rhoinf, pinf);
                    double k = spec_kine_energy(vinf1(0), vinf1(1), vinf1(2));
                    double E = total_energy(rhoinf, k, e);
                    double H = total_enthalpy(rhoinf, pinf, E);
                    double a = GAMMA + 1.;
                    double b = (6.-2.*GAMMA)*rminus;
                    double cc = (GAMMA + 1.) * std::pow(rminus, 2.) - 16. * H;

                    double rplus1 = (-b - std::sqrt(std::pow(b, 2.) - 4. * a * cc)) / (2. * a);
                    double rplus2 = (-b + std::sqrt(std::pow(b, 2.) - 4. * a * cc)) / (2. * a);
                    if (std::pow(b, 2.) - 4. * a * cc < 0.)
                    {
                        mc.prim_[0] = rhoinf;
                        mc.prim_[1] = vinf1(0);
                        mc.prim_[2] = vinf1(1);
                        mc.prim_[3] = vinf1(2);
                        mc.prim_[4] = pinf;

                        mc.cons_sp1_ = prim_to_cons(mc.prim_);

                        assert(p >= 0.);
                        assert(rho > 0.);

                        continue;
                    }
                    //std::cout << "pbound: " << pbound << std::endl;
                    //std::cout << "rhobound: " << rhobound << std::endl;
                    //std::cout << "vbound: " << vbound << std::endl;
                    //std::cout << "cbound: " << cbound << std::endl;
                    //std::cout << "rminus: " << rminus << std::endl;
                    //std::cout << "rplus1: " << rplus1 << std::endl;
                    //std::cout << "rplus2: " << rplus2 << std::endl;

                    rplus = rplus1;
                    v = 0.5 * (rplus + rminus);
                    Vector3 vv;
                    vv = normal * v;
                    double c = 0.25 * (rplus - rminus) * (GAMMA - 1.);

                    s = sinf;
                    rho = std::pow(s*c*c/GAMMA, 1. / (GAMMA - 1.));
                    p = rho * c * c / GAMMA;
                }
                else
                {
                    assert(false);
                    subsonic = false;
                    mc.prim_[0] = rhoinf;
                    mc.prim_[1] = vinf1(0);
                    mc.prim_[2] = vinf1(1);
                    mc.prim_[3] = vinf1(2);
                    mc.prim_[4] = pinf;

                    mc.cons_sp1_ = prim_to_cons(mc.prim_);

                    continue;
                }
            }
            else // outflow
            {
                inflow = false;
                if (dotp(vinf, normal) < cinf) // subsonic
                {
                    subsonic = true;
                    //std::cout << "subsonic outflow" << std::endl;
                    double c = cbound;
                    rplus = vbound + 2. * cbound / (GAMMA - 1.);
                    rminus = rplus - 4. * cbound / (GAMMA - 1.);
                    s = sbound;

                    v = 0.5 * (rplus + rminus);

                    rho = std::pow(s*c*c/GAMMA, 1. / (GAMMA - 1.));
                    p = rho * c * c / GAMMA;
                }
                else
                {
                    assert(false);
                    subsonic = false;
                    Vector3 vbound2 = normal * dotp(vbound1, normal);
                    mc.prim_[0] = rhobound;
                    mc.prim_[1] = vbound2(0);
                    mc.prim_[2] = vbound2(1);
                    mc.prim_[3] = vbound2(2);
                    mc.prim_[4] = pbound;

                    mc.cons_sp1_ = prim_to_cons(mc.prim_);

                    continue;
                }
            }

            //double v = 0.5 * (rplus + rminus);
            //double c = 0.25 * (rplus - rminus) * (GAMMA - 1.);
            //double s;
            Vector3 vv;
            vv = normal * v;

            //if (vbound > 0.) // outflow
            //{
                //vv = vbound1 + normal * (v - vbound); 
                //s = sbound;
            //}
            //else
            //{
                //vv = vinf1 + normal * (v - vinf.len());
                //s = sinf;
            //}

            //double rho = std::pow(s*c*c/GAMMA, 1. / (GAMMA - 1.));
            //double p = rho * c * c / GAMMA;

            mc.prim_[0] = rho;
            mc.prim_[1] = vv(0);
            mc.prim_[2] = vv(1);
            mc.prim_[3] = vv(2);
            mc.prim_[4] = p;

            if (subsonic)
            {
                //if (std::abs(vv(0)) > std::abs(vinf(0)))
                //{
                    //std::cout << "vv(0): " << vv(0) << std::endl;
                //}
                //assert(vv.len() - vinf1.len() <= ZERO);
            }

            //if (inflow)
            {
            //std::cout << "inflow: " << inflow << std::endl;
            //std::cout << "rho: " << rho << std::endl;
            //std::cout << "u: " << vv(0) << std::endl;
            //std::cout << "v: " << vv(1) << std::endl;
            //std::cout << "w: " << vv(2) << std::endl;
            //std::cout << "p: " << p << std::endl;
            }

            mc.cons_sp1_ = prim_to_cons(mc.prim_);

            assert(rho > 0.);
            assert(p >= 0.);
        }
    }*/

    /*void BoundaryCondition::init_sod(Mesh& mesh)
    {
        std::ifstream in;
        in.open("init_sod.bc");
        if (!in.is_open()) {
            return;
        }

        Vector5 primL, primR;
        double centerx;

        in >> primL[0];
        in >> primL[1];
        in >> primL[2];
        in >> primL[3];
        in >> primL[4];

        in >> primR[0];
        in >> primR[1];
        in >> primR[2];
        in >> primR[3];
        in >> primR[4];
        
        in >> centerx;

        in.close();

        for (MeshCell& mc: mesh.cell_)
        {
            if (mc.poly().centroid()(0) <= centerx)
            {
                mc.prim_ = primL;
            }
            else
            {
                mc.prim_ = primR;
            }

            mc.cons_sp1_ = prim_to_cons(mc.prim_);
        }

        for (MeshCell& mc: mesh.dirichlet_boundaries_)
        {
            if (mc.poly().centroid()(0) <= centerx)
            {
                mc.prim_ = primL;
            }
            else
            {
                mc.prim_ = primR;
            }

            mc.cons_sp1_ = prim_to_cons(mc.prim_);
        }
    }*/

    /*void BoundaryCondition::farfield(Mesh& mesh)
    {
        // https://github.com/su2code/SU2_EDU/blob/master/src/solver_direct_mean.cpp#L2384

        std::ifstream in;
        in.open("init_farfield.bc");
        if (!in.is_open()) {
            return;
        }

        double rhoinf;
        double pinf;
        double machinf;
        double aoa;

        in >> rhoinf;
        in >> pinf;
        in >> machinf;
        in >> aoa;

        in.close();

        //double vb = 0.;
        //double rhoinf = 1.;
        //double pinf = 0.71428571428;
        //double machinf = 0.8;
        //double aoa = 1.25;

        double cinf = std::sqrt(GAMMA * pinf / rhoinf);
        //double vinf = machinf * cinf;
        Vector3 vinf(machinf * cinf * std::cos(deg_to_rad(aoa)), machinf * cinf * std::sin(deg_to_rad(aoa)), 0.);
        double sinf = std::pow(rhoinf, GAMMA) / pinf;

        for (MeshCell& mc: mesh.farfield_boundaries_)
        {
            assert(mc.face().size() == 1);
            Vector3 normal = mc.face()[0].face().normal();

            //Vector3 vinf1 = normal * vinf;
            Vector3 vinf1 = normal * dotp(vinf, normal);

            const MeshCell& interior = mesh.cell(mc.interior_boundary());

            double rhobound = interior.prim(0);
            double pbound = interior.prim(4);
            Vector3 vbound1(interior.prim(1), interior.prim(2), interior.prim(3));
            double vbound = dotp(vbound1, normal);
            double cbound = std::sqrt(GAMMA * pbound / rhobound);
            double sbound = std::pow(rhobound, GAMMA) / pbound;

            double rplus, rminus;
            if (dotp(vinf, normal) > -cinf)
            //if (vinf1.len() > -cinf)
            {
                //--- Subsonic inflow or outflow ---
                rplus = vbound + 2. * cbound / (GAMMA - 1.);
            }
            else
            {
                //--- Supersonic inflow ---
                //rplus = vinf1.len() + 2. * cinf / (GAMMA - 1.);
                rplus = dotp(vinf, normal) + 2. * cinf / (GAMMA - 1.);
            }

            if (dotp(vinf, normal) > cinf)
            //if (vinf1.len() > cinf)
            {
                //--- Supersonic outflow ---
                rminus = vbound - 2. * cbound / (GAMMA - 1.);
            }
            else
            {
                //--- Subsonic outflow ---
                //rminus = vinf.len() - 2. * cinf / (GAMMA - 1.);
                rminus = dotp(vinf, normal) - 2. * cinf / (GAMMA - 1.);
            }

            double v = 0.5 * (rplus + rminus);
            double c = 0.25 * (rplus - rminus) * (GAMMA - 1.);
            double s;
            Vector3 vv;


            if (dotp(vinf, normal) > 0.)
            //if (vinf1.len() > 0.)
            {
                //--- Outflow conditions ---
                vv = vbound1 + normal * (v - vbound); 
                s = sbound;
            }
            else
            {
                //--- Inflow conditions ---
                vv = vinf1 + normal * (v - vinf.len());
                s = sinf;
            }

            double rho = std::pow(s*c*c/GAMMA, 1. / (GAMMA - 1.));
            double p = rho * c * c / GAMMA;

            mc.prim_[0] = rho;
            mc.prim_[1] = vv(0);
            mc.prim_[2] = vv(1);
            mc.prim_[3] = vv(2);
            mc.prim_[4] = p;

            //std::cout << "rho: " << rho << std::endl;
            //std::cout << "u: " << vv(0) << std::endl;
            //std::cout << "v: " << vv(1) << std::endl;
            //std::cout << "w: " << vv(2) << std::endl;
            //std::cout << "p: " << p << std::endl;

            mc.cons_sp1_ = prim_to_cons(mc.prim_);

            assert(rho > 0.);
            assert(p > 0.);
        }
    }*/

    void BoundaryCondition::empty(Mesh& mesh, MeshCell& mc)
    {
        assert(mesh.query(mc.interior_boundary()) != nullptr); 
        mc.prim_ = mesh.cell(mc.interior_boundary()).prim();
        mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);
    }

    void BoundaryCondition::symmetry(Mesh& mesh, MeshCell& mc)
    {
        // Implement as slip wall.
        slipwall(mesh, mc);
    }

    void BoundaryCondition::slipwall(Mesh& mesh, MeshCell& mc)
    {
        const MeshCell& nei = mesh.cell(mc.interior_boundary());

        //mc.prim_[0] = nei.prim(0);
        //mc.prim_[1] = 0.;
        //mc.prim_[2] = nei.prim(2);
        //mc.prim_[3] = nei.prim(3);
        //mc.prim_[4] = nei.prim(4);

        //assert(mc.prim_[0] > 0.);
        //mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);

        //return;

        //const MeshFace& mf = mesh.face(mc.face()[0].tag());
        auto mf = std::find_if(nei.face().begin(), nei.face().end(), [&](const MeshFace& f){return f.tag() == mc.face()[0].tag();});
        if (mf == nei.face().end())
        {
            std::cout << "mesh: " << mesh.tag()() << std::endl;
            std::cout << "mf: " << mc.face()[0].tag() << std::endl;
            //for (const auto& f: nei.face())
            for (auto f = nei.face().begin(); f != nei.face().end(); ++f)
            {
                std::cout << "nei mf: " << f->tag() << std::endl;
            }
        }
        assert(mf != nei.face().end());
        //const MeshFace& mf = mc.face()[0];
        auto n = mf->face().normal();

        //if (mesh.tag()() == 1 && nei.tag()() == 8579)
        //if (mesh.tag()() == 1)
        //{
            //std::cout << "n(0): " << n(0) << std::endl;
            //std::cout << "n(1): " << n(1) << std::endl;
            //std::cout << "n(2): " << n(2) << std::endl;
            //assert(false);
        //}
        //for (const MeshFace& f: nei.face())
        //{
        //    if (f.tag() == mf->tag())
        //    {
        //        assert(dotp(f.face().normal(), n) > TAILOR_ZERO);
        //        break;
        //    }
        //}
        assert(dot(mf->face().normal(), n) > TAILOR_ZERO);

        auto T = make_rot_matrix(n);
        auto neii = nei.prim();
        auto vf = mf->vf();
        neii(1) -= vf(0);
        neii(2) -= vf(1);
        neii(3) -= vf(2);
        //auto prim = T * nei.prim();
        auto prim = T * neii;
        prim(1) *= -1.;

        mc.prim_ = T.transpose() * prim;
        mc.prim_(1) += vf(0);
        mc.prim_(2) += vf(1);
        mc.prim_(3) += vf(2);
        mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);

        return;

        mc.prim_(0) = nei.prim(0);
        mc.prim_(4) = nei.prim(4);

        //if (mc.tag()() == 40)
        //{
            //std::cout << "aa: " << mc.prim
            //assert(false);
        //}
        //assert(mc.prim_[0] < 1e9);

        Vector3 vel(nei.prim(1), nei.prim(2), nei.prim(3));

        assert(!nei.prim_.isnan());

        auto normvel = n * dot(vel, n);
        //auto tangvel = vel - normvel;

        //assert(!std::isnan(mf->vgn()));

        //double vgx = mf->vgn() * n(0);
        //double vgy = mf->vgn() * n(1);
        //double vgz = mf->vgn() * n(2);
        //double vgx = mf->vgn()(0);
        //double vgy = mf->vgn()(1);
        //double vgz = mf->vgn()(2);
        auto aa = n * mf->vgn();
        double vgx = aa(0);
        double vgy = aa(1);
        double vgz = aa(2);
        Vector3 vgn(vgx, vgy, vgz);

        //mc.prim_[1] = tangvel(0) + vgx;
        //mc.prim_[2] = tangvel(1) + vgy;
        //mc.prim_[3] = tangvel(2) + vgz;

        Vector3 v = vel - normvel * 2. + vgn;
        mc.prim_(1) = v(0);
        mc.prim_(2) = v(1);
        mc.prim_(3) = v(2);

        mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);
    }

    void BoundaryCondition::init_farfield(Mesh& mesh, const Vector5& prim)
    {
        for (MeshCell& mc: mesh.farfield_boundaries_)
        {
            mc.prim_ = prim;
            mc.cons_sp1_ = prim_to_cons(prim, fs_.gamma_);
        }
    }

    void BoundaryCondition::set_dirichlet(Mesh& mesh)
    {
        if (mesh.dirichlet_boundaries_.empty())
        {
            return;
        }

        //auto prim = read_dirichlet(mesh.tag());
        // prim should already be set by init_flow.

        for (MeshCell& mc: mesh.dirichlet_boundaries_)
        {
            const MeshCell& nei = mesh.cell(mc.interior_boundary());
            auto mf = std::find_if(nei.face().begin(), nei.face().end(), [&](const MeshFace& f){return f.tag() == mc.face()[0].tag();});

            // TODO WHy add face vel???
            //mc.prim_ = prim;
            auto vf = mf->vf();
            mc.prim_(1) -= vf(0);
            mc.prim_(2) -= vf(1);
            mc.prim_(3) -= vf(2);
            mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);
        }
    }

    //void BoundaryCondition::set_dirichlet(Mesh& mesh, const Vector5& prim)
    //{
        //for (MeshCell& mc: mesh.dirichlet_boundaries_)
        //{
            //mc.prim_ = prim;
            //mc.cons_sp1_ = prim_to_cons(prim, fs_.gamma_);
        //}
    //}

    void BoundaryCondition::interog(Mesh& mesh)
    {
        for (MeshCell& mc: mesh.interog_boundaries_)
        {   
            const MeshCell& nei = mesh.cell(mc.interior_boundary());
            const auto& mf = mc.face()[0];
            //const auto& vgn = mf.vgn();
            const auto& vgn = mf.face().normal() * mf.vgn();

            mc.prim_ = nei.prim();

            mc.prim_(1) -= vgn(0);
            mc.prim_(2) -= vgn(1);
            mc.prim_(3) -= vgn(2);

            mc.cons_sp1_ = prim_to_cons(mc.prim_, fs_.gamma_);
        }
    }

    void BoundaryCondition::update_fs(const Freestream& fs)
    {
        fs_ = fs;
    }

    void BoundaryCondition::set_bc(Mesh& mesh, Profiler* profiler)
    {       
        //if (profiler != nullptr) {profiler->start("bc-read");}
        if (first_)
        {
            fs_.read();
            first_ = false;

            set_dirichlet(mesh);
        }
        //if (profiler != nullptr) {profiler->stop("bc-read");}

        interog(mesh);

        //if (profiler != nullptr) {profiler->start("bc-wall");}
        for (MeshCell& mc: mesh.wall_boundaries_)
        {   
            slipwall(mesh, mc);
            //assert(mc.prim(0) > 0.);
        }
        for (MeshCell& mc: mesh.symmetry_boundaries_)
        {   
            symmetry(mesh, mc);
        }
        //if (profiler != nullptr) {profiler->stop("bc-wall");}

        //if (profiler != nullptr) {profiler->start("bc-empty");}
        for (MeshCell& mc: mesh.empty_boundaries_)
        {   
            empty(mesh, mc);
            //assert(mc.prim(0) > 0.);
        }
        //if (profiler != nullptr) {profiler->stop("bc-empty");}

        //if (profiler != nullptr) {profiler->start("bc-far");}
        //farfield5(mesh);
        farfield4(mesh); // TODO try this
        //if (profiler != nullptr) {profiler->stop("bc-far");}
    }

    Vector5 BoundaryCondition::read_dirichlet(const Tag& meshtag)
    {
        namespace po = boost::program_options;
        po::options_description op;

        std::string sdesc = "mesh ";
        sdesc.append(std::to_string(meshtag()));

        po::options_description desc{sdesc};
        desc.add_options()
            (cstr(sdesc, "rho"), po::value<double>()->required(), "")
            (cstr(sdesc, "p"), po::value<double>()->required(), "")
            (cstr(sdesc, "u"), po::value<double>()->required(), "")
            (cstr(sdesc, "v"), po::value<double>()->required(), "")
            (cstr(sdesc, "w"), po::value<double>()->required(), "")
            ;

        op.add(desc);
        std::string fn = "dirichlet.ini";
        std::ifstream settings_file(fn);

        boost::program_options::variables_map vm;
        po::store(po::parse_config_file(settings_file, op, true), vm);
        po::notify(vm);

        double rho = vm[cstr(sdesc, "rho")].as<double>();
        double p = vm[cstr(sdesc, "p")].as<double>();
        double u = vm[cstr(sdesc, "u")].as<double>();
        double v = vm[cstr(sdesc, "v")].as<double>();
        double w = vm[cstr(sdesc, "w")].as<double>();

        Vector5 prim;
        prim(0) = rho;
        prim(1) = u;
        prim(2) = v;
        prim(3) = w;
        prim(4) = p;

        return prim;
    }
}
