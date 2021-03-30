#include "riemann.h"

namespace Tailor
{
    RiemannSolver::RiemannSolver(RiemannSolverType riemann_solver_type, const State& left_state, const State& right_state, double face_area, double gamma, double& max_eigen, Vector5& flux, SpeedEstimateHLLC sehllc): speed_estimate_hllc_(sehllc)
    {
        if (riemann_solver_type == RiemannSolverType::roe)
        {
            Matrix5 Aroe;
            roe(left_state, right_state, flux, Aroe, max_eigen, face_area, gamma);

            assert(!flux.isnan());
            assert(max_eigen > 0.);
        }
        else if (riemann_solver_type == RiemannSolverType::hllc)
        {
            double sL, sR;
            std::tie(flux, max_eigen, sL, sR) = hllc(left_state, right_state, face_area, gamma);
        }
        else
        {
            assert(false);
        }
    }
    
    Vector5 RiemannSolver::jump_ver_2(const State& left, const State& right)
    {
        double rhoL = left.rho;
        double uL = left.u;
        double vL = left.v;
        double wL = left.w;
        double EL = left.E;
        double rhoR = right.rho;
        double uR = right.u;
        double vR = right.v;
        double wR = right.w;
        double ER = right.E;

        Vector5 jump;

        jump(0) = rhoR - rhoL;
        jump(1) = rhoR*uR - rhoL*uL; 
        jump(2) = rhoR*vR - rhoL*vL; 
        jump(3) = rhoR*wR - rhoL*wL; 
        jump(4) = ER - EL; 

        return jump;
    }

    //Vector5 RiemannSolver::jump_ver_3(const State& left, const State& right)
    //{
    //    double rhoL = left.rho;
    //    double uL = left.u;
    //    double vL = left.v;
    //    double wL = left.w;
    //    double pL = left.p;
    //    double rhoR = right.rho;
    //    double uR = right.u;
    //    double vR = right.v;
    //    double wR = right.w;
    //    double pR = right.p;

    //    Vector5 jump;

    //    jump[0] = rhoR - rhoL;
    //    jump[1] = uR - uL;
    //    jump[2] = vR - vL;
    //    jump[3] = wR - wL;
    //    jump[4] = pR - pL;

    //    return jump;
    //}

    void RiemannSolver::calc_roe_ave_vars(const State& left, const State& right, double gamma)
    {
        double rhoL = left.rho;
        double uL = left.u;
        double vL = left.v;
        double wL = left.w;
        double HL = left.H;
        double rhoR = right.rho;
        double uR = right.u;
        double vR = right.v;
        double wR = right.w;
        double HR = right.H;

        double srL = std::sqrt(rhoL);
        double srR = std::sqrt(rhoR);
        double denom = srL + srR;
        u = (srL * uL + srR * uR) / denom;
        v = (srL * vL + srR * vR) / denom;
        w = (srL * wL + srR * wR) / denom;
        H = (srL * HL + srR * HR) / denom;

        //double RT = std::sqrt(rhoR/rhoL);
        //rho = RT * rhoL;
        //u = (uL + RT * uR) / (1. + RT);
        //v = (vL + RT * vR) / (1. + RT);
        //w = (wL + RT * wR) / (1. + RT);
        //H = (HL + RT * HR) / (1. + RT);

        k = spec_kine_energy(u, v, w);
        a = std::sqrt((gamma - 1.) * (H - k));

        if (H <= k || std::isnan(H) || std::isnan(k))
        {
            std::cout << "LEFT VALUES" << std::endl;
            left.print();
            std::cout << "RIGHT VALUES" << std::endl;
            right.print();
            std::cout << "ROE VALUES" << std::endl;
            std::cout << "H: " << H << std::endl;
            std::cout << "k: " << k << std::endl;
        }
        assert(H > k);
        assert(!std::isnan(H));
        assert(!std::isnan(k));
        assert(!std::isnan(a));
        assert(!std::isinf(a));

        //double nx = left.nx;
        //double ny = left.ny;
        //double nz = left.nz;

        //double lx = left.lx;
        //double ly = left.ly;
        //double lz = left.lz;

        //double mx = left.mx;
        //double my = left.my;
        //double mz = left.mz;

        //qn  = u*nt.n(0) + v*nt.n(1) + w*nt.n(2);
        //ql  = u*nt.l(0) + v*nt.l(1) + w*nt.l(2);
        //qm  = u*nt.m(0) + v*nt.m(1) + w*nt.m(2);
        //qn  = u*nx + v*ny + w*nz;
        //ql  = u*lx + v*ly + w*lz;
        //qm  = u*mx + v*my + w*mz;
        //asq = std::pow(a,2.);
    }

    Matrix5 RiemannSolver::right_eigenv()
    {
        Matrix5 R;

        double vnet = u;

        R(0,0) = 1.;
        R(1,0) = vnet - a;
        R(2,0) = v;
        R(3,0) = w;
        R(4,0) = H - vnet * a;

        R(0,1) = 1.;
        R(1,1) = vnet;
        R(2,1) = v;
        R(3,1) = w;
        R(4,1) = k;

        R(0,2) = 0.;
        R(1,2) = 0.;
        R(2,2) = 1.;
        R(3,2) = 0.;
        R(4,2) = v;

        R(0,3) = 0.;
        R(1,3) = 0.;
        R(2,3) = 0.;
        R(3,3) = 1.;
        R(4,3) = w;

        R(0,4) = 1.;
        R(1,4) = vnet + a;
        R(2,4) = v;
        R(3,4) = w;
        R(4,4) = H + vnet * a;

        return R;
    }

    void rhll_normals(const State& left, const State& right, const Vector3& normal, double& alpha1, double& alpha2, Vector3& n1, Vector3& n2)
    {
        double nx = normal(0);
        double ny = normal(1);
        double nz = normal(2);

        double uL = left.u;
        double vL = left.v;
        double wL = left.w;

        double uR = right.u;
        double vR = right.v;
        double wR = right.w;

        double abs_dq = std::sqrt(std::pow((uR-uL), 2.) + std::pow((vR-vL), 2.) + std::pow((wR-wL), 2.));

        double nx1, ny1, nz1;
        double ss;
        double xx = -1;
        double yy = -1;
        double zz = -1;
        if (abs_dq > TAILOR_ZERO)
        {
            nx1 = (uR-uL) / abs_dq;
            ny1 = (vR-vL) / abs_dq;
            nz1 = (wR-wL) / abs_dq;
        }
        else
        {
            double tempx = ny*ny + nz*nz;
            double tempy = nz*nz + nx*nx;
            double tempz = nx*nx + ny*ny;

            if (tempx >= tempy && tempx >= tempz)
            {
                nx1 = 0.;
                ny1 = -nz;
                nz1 = ny;
            }
            else if (tempy >= tempx && tempy >= tempz)
            {
                nx1 = -nz;
                ny1 = 0.;
                nz1 = nx;
            }
            else if (tempz >= tempx && tempz >= tempy)
            {
                nx1 = -ny;
                ny1 = nx;
                nz1 = 0.;
            }
            else
            {
                assert(false);
            }
            
            xx = nx1;
            yy = ny1;
            zz = nz1;

            if (xx == 0. && yy == 0. && zz == 0.)
            {
                std::cout << "nx: " << nx << std::endl;
                std::cout << "ny: " << ny << std::endl;
                std::cout << "nz: " << nz << std::endl;
                std::cout << "tempx: " << tempx << std::endl;
                std::cout << "tempy: " << tempy << std::endl;
                std::cout << "tempz: " << tempz << std::endl;
                assert(false);
            }

            double temp = std::sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);
            nx1 = nx1 / temp;
            ny1 = ny1 / temp;
            nz1 = nz1 / temp;

        }

        alpha1 = nx*nx1 + ny*ny1 + nz*nz1;

        ss = sign(1., alpha1);
        double temp = sign(1., alpha1);
        nx1 = temp * nx1;
        ny1 = temp * ny1;
        nz1 = temp * nz1;
        alpha1 = temp * alpha1;

        n1(0) = nx1;
        n1(1) = ny1;
        n1(2) = nz1;

        double tempx = ny1*nz - nz1*ny;
        double tempy = nz1*nx - nx1*nz;
        double tempz = nx1*ny - ny1*nx;

        double nx2 = tempy*nz1 - tempz*ny1;
        double ny2 = tempz*nx1 - tempx*nz1;
        double nz2 = tempx*ny1 - tempy*nx1;

        temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 );
        if (temp == 0.)
        {
            std::cout << "nx1: " << nx1 << std::endl;
            std::cout << "ny1: " << ny1 << std::endl;
            std::cout << "nz1: " << nz1 << std::endl;
            std::cout << "nx2: " << nx2 << std::endl;
            std::cout << "ny2: " << ny2 << std::endl;
            std::cout << "nz2: " << nz2 << std::endl;
            std::cout << "nx: " << nx << std::endl;
            std::cout << "ny: " << ny << std::endl;
            std::cout << "nz: " << nz << std::endl;
            std::cout << "tempx: " << tempx << std::endl;
            std::cout << "tempy: " << tempy << std::endl;
            std::cout << "tempz: " << tempz << std::endl;
            std::cout << "alpha1: " << alpha1 << std::endl;
            std::cout << "abs_dq: " << abs_dq << std::endl;
            std::cout << "ss: : " << ss << std::endl;
            std::cout << "xx: : " << xx << std::endl;
            std::cout << "yy: : " << yy << std::endl;
            std::cout << "zz: : " << zz << std::endl;
        }
        assert(temp != 0.);
        nx2 = nx2 / temp;
        ny2 = ny2 / temp;
        nz2 = nz2 / temp;

        alpha2 = nx*nx2 + ny*ny2 + nz*nz2;

        temp = sign(1. , alpha2);
        nx2 = temp * nx2;
        ny2 = temp * ny2;
        nz2 = temp * nz2;
        alpha2 = temp * alpha2;

        n2(0) = nx2;
        n2(1) = ny2;
        n2(2) = nz2;

        assert(!n1.isnan());
        assert(!n2.isnan());
    }

    std::tuple<double, double> RiemannSolver::rhll_ws_est_pres(const State& left, const State& right, double gamma)
    {
        double rhoL = left.rho;
        double uL = left.u;
        double vL = left.v;
        double wL = left.w;
        double aL = left.a;
        double pL = left.p;

        double rhoR = right.rho;
        double uR = right.u;
        double vR = right.v;
        double wR = right.w;
        double aR = right.a;
        double pR = right.p;

        double rhodash = 0.5 * (rhoL + rhoR);
        double adash   = 0.5 * (aL + aR);
        double ppvrs = 0.5 * (pL + pR) - 0.5 * (uR - uL) * rhodash * adash;
        double ps = std::max(0., ppvrs);

        double qL, qR;

        if (ps <= pL)
        {
            qL = 1.;
        }
        else
        {
            qL = std::sqrt(1. + (gamma + 1.) * (ps / pL - 1.) / (2. * gamma));
        }

        if (ps <= pR)
        {
            qR = 1.;
        }
        else
        {
            qR = std::sqrt(1. + (gamma + 1.) * (ps / pR - 1.) / (2. * gamma));
        }

        double sL = uL - aL * qL;
        double sR = uR + aR * qR;

        return std::make_tuple(sL, sR);
    }

    std::tuple<double, double> RiemannSolver::rhll_ws_est(const State& left, const State& right)
    {
        double uL = left.u;
        double vL = left.v;
        double wL = left.w;
        double aL = left.a;

        double uR = right.u;
        double vR = right.v;
        double wR = right.w;
        double aR = right.a;

        //SLm = std::min({0., u - a + vfn, uL - aL + vfn});
        //SRp = std::max({0., u + a + vfn, uR + aR + vfn});

        double sL = std::min({0., u - a, uL - aL});
        double sR = std::max({0., u + a, uR + aR});

        //SLm = std::min({u - a, uL - aL});
        //SRp = std::max({u + a, uR + aR});
        
        //SLm = u - a;
        //SRp = u + a;

        return std::make_tuple(sL, sR);
    }

    /*Matrix5 RiemannSolver::rhll_abs_eigen(double alpha1, double alpha2, double SLm, double SRp, double vfn)
    {
        auto ws = abs_eigen(vfn);
        auto eig = eigen(vfn);
        Matrix5 I;
        I = 0.;
        double temp = alpha1 * 2. * SRp * SLm;
        I(0,0) = temp;
        I(1,1) = temp;
        I(2,2) = temp;
        I(3,3) = temp;
        I(4,4) = temp;

        auto comws = ws * alpha2 - (I + eig * alpha2 * (SRp + SLm)) * (1. / (SRp - SLm));

        return comws;
    }*/

    Matrix5 RiemannSolver::eigen()
    {
        // wave speeds

        Matrix5 ws;
        ws = 0.;

        assert(!std::isnan(a));

        double vnet = u;

        ws(0,0) = vnet - a;
        ws(1,1) = vnet;
        ws(2,2) = vnet;
        ws(3,3) = vnet;
        ws(4,4) = vnet + a;

        return ws;
    }

    Matrix5 RiemannSolver::abs_eigen(const State& left, const State& right)
    {
        // Absolute values of the wave speeds

        Matrix5 ws;
        ws = 0.;

        assert(!std::isnan(a));

        double vnet = u;

        ws(0,0) = std::abs(vnet - a);
        ws(1,1) = std::abs(vnet);
        ws(2,2) = std::abs(vnet);
        ws(3,3) = std::abs(vnet);
        ws(4,4) = std::abs(vnet + a);

        //double lL[5];
        //double lR[5];
        //double l[5];

        //lL[0] = left.u - left.a;
        //lL[1] = left.u;
        //lL[2] = left.u;
        //lL[3] = left.u;
        //lL[4] = left.u + left.a;

        //lR[0] = right.u - right.a;
        //lR[1] = right.u;
        //lR[2] = right.u;
        //lR[3] = right.u;
        //lR[4] = right.u + right.a;

        //l[0] = u - a;
        //l[1] = u;
        //l[2] = u;
        //l[3] = u;
        //l[4] = u + a;

        //for (int i=0; i<5; ++i)
        //{
        //    double sigma = std::max(0., l[i] - lL[i]);
        //    sigma = std::max(sigma, lR[i] - l[i]);

        //    double akL = l[i] - sigma; 
        //    double akR = l[i] + sigma; 

        //    //if (lL[i] < 0. && lR[i] > 0.)
        //    if (akL < 0. && akR > 0.)
        //    {
        //        //ws(i,i) = ((lR[i] + lL[i]) * l[i] - 2. * lR[i] * lL[i]) / (lR[i] - lL[i]);
        //        ws(i,i) = ((akR + akL) * l[i] - 2. * akR * akL) / (akR - akL);
        //    }
        //}

        return ws;
    }

    /*Vector5 RiemannSolver::wave_strength_ver_2(const State& left, const State& right, const Vector3& n)
    {
        auto jump = jump_ver_2(left, right);

        Vector5 ldu;

        double nx = n(0);
        double ny = n(1);
        double nz = n(2);

        double B = jump[2] - v * jump[0];
        double A = (jump[3] - w * jump[0] - nz * B / ny) / nx;

        ldu[3] = ny * (jump[1] - u * jump[0] - nx * B / ny + nz * A) / (std::pow(nx, 2.) + std::pow(ny, 2.) + std::pow(nz, 2.));
        ldu[4] = A - nz * ldu[3] / ny;
        ldu[1] = (jump[4] - H * jump[0] - qn * B / ny - ldu[3] * (qn * nx / ny + u * ny - v * nx) - ldu[4] * (w * nx - u * nz)) / (k - H);
        ldu[2] = 0.5 * (jump[0] - ldu[1] + (B + ldu[3] * nx) / (a * ny));
        ldu[0] = jump[0] - ldu[1] - ldu[2];

        return ldu;
    }*/

    Vector5 RiemannSolver::wave_strength_ver_1(const State& left, const State& right, double gamma)
    {
        Vector5 jump = jump_ver_2(left, right);
        Vector5 ldu;
        
        auto vnet = u;

        ldu(2) = jump(2) - v * jump(0);
        ldu(3) = jump(3) - w * jump(0);

        double j4 = jump(4) - ldu(2) * v - ldu(3) * w;

        ldu(1) = (gamma - 1.) * (jump(0) * (H - std::pow(vnet, 2.)) + vnet * jump(1) - j4) / std::pow(a, 2.); 
        ldu(0) = (jump(0) * (vnet + a) - jump(1) - a * ldu(1)) / (2. * a); 
        ldu(4) = jump(0) - (ldu(0) + ldu(1));

        assert(!jump.isnan());
        

        //double rhoL = left.rho;
        //double vL = left.u;
        //double pL = left.p;

        //double rhoR = right.rho;
        //double vR = right.u;
        //double pR = right.p;

        //double drho = rhoR - rhoL;
        //double du =   vR - vL;
        //double dP =   pR - pL;

        //double RT = sqrt(rhoR/rhoL);
        //double rho = RT*rhoL;

        //ldu[0] =  0.5*(dP-rho*a*du)/(a*a);
        //ldu[1] = -( dP/(a*a) - drho );
        //ldu[4] =  0.5*(dP+rho*a*du)/(a*a);

        return ldu;
    }

    Vector5 RiemannSolver::dissipation_term(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double gamma)
    {
        Vector5 diss;

        Vector5 ldu = wave_strength_ver_1(left, right, gamma);
        //auto ldu = wave_strength_ver_2(left, right, n);


        for (int i=0; i<NVAR; ++i)
        {
            diss(i) = 0.;
            for (int j=0; j<NVAR; ++j)
            {
                assert(!std::isnan(ws(j,j)));
                assert(!std::isnan(ldu(j)));
                assert(!std::isnan(R(i,j)));
                diss(i) += ws(j,j) * ldu(j) * R(i,j);
            }
        }

        return diss;
    }

    Vector5 RiemannSolver::rhll_numerical_flux(double SLm, double SRp, const State& left, const State& right, const State& leftn2, const State& rightn2, const Matrix5& ws, const Matrix5& R, double facearea, const Vector3& n, double gamma)
    {
        Vector5 diss = dissipation_term(leftn2, rightn2, ws, R, gamma);
        auto nf = (SRp*left.flux - SLm*right.flux)/(SRp-SLm) - 0.5 * diss;
        nf = facearea * nf;
        return nf;
    }

    Vector5 RiemannSolver::hlle_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea)
    {
        auto nf = (SRp*left.flux - SLm*right.flux + SLm*SRp*(right.cons-left.cons))/(SRp-SLm);
        nf = facearea * nf;
        return nf;
    }

    Vector5 RiemannSolver::hllc_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea)
    {
        auto rhoL = left.rho;
        auto uL = left.u;
        auto vL = left.v;
        auto wL = left.w;
        auto pL = left.p;
        auto EL = left.E;
        auto fL = left.flux;

        auto rhoR = right.rho;
        auto uR = right.u;
        auto vR = right.v;
        auto wR = right.w;
        auto pR = right.p;
        auto ER = right.E;
        auto fR = right.flux;

        auto temp1 = rhoL * (SLm - uL);
        auto temp2 = rhoR * (SRp - uR);

        auto Ss = (pR - pL + uL * temp1 - uR * temp2) / (temp1 - temp2);

        Vector5 vecL;
        vecL(0) = 1.;
        vecL(1) = Ss;
        vecL(2) = vL;
        vecL(3) = wL;
        vecL(4) = EL / rhoL + (Ss - uL) * (Ss + pL / (rhoL * (SLm - uL)));

        Vector5 vecR;
        vecR(0) = 1.;
        vecR(1) = Ss;
        vecR(2) = vR;
        vecR(3) = wR;
        vecR(4) = ER / rhoR + (Ss - uR) * (Ss + pR / (rhoR * (SRp - uR)));

        auto usL = vecL * rhoL * ((SLm - uL) / (SLm - Ss));
        auto usR = vecR * rhoR * ((SRp - uR) / (SRp - Ss));

        auto fsL = fL + (usL - left.cons) * SLm;
        auto fsR = fR + (usR - right.cons) * SRp;

        if (SLm >= 0.)
        {
            return facearea * fL;
        }
        else if (SLm <= 0. && Ss >= 0.)
        {
            return facearea * fsL;
        }
        else if (Ss <= 0. && SRp >= 0.)
        {
            return facearea * fsR;
        }
        else if (SRp <= 0.)
        {
            return facearea * fR;
        }
        else
        {
            std::cout << "SL: " << SLm << std::endl;
            std::cout << "SR: " << SRp << std::endl;
            std::cout << "Ss: " << Ss << std::endl;
            assert(false);
        }

        //if (Ss >= 0.)
        //{
        //    return facearea * fsL;
        //}
        //else
        //{
        //    return facearea * fsR;
        //}
    }

    Vector5 RiemannSolver::numerical_flux(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double facearea, double gamma)
    {
        Vector5 diss = dissipation_term(left, right, ws, R, gamma);
        assert(!left.flux.isnan());
        assert(!right.flux.isnan());
        assert(!diss.isnan());
        assert(!std::isnan(facearea));
        //auto nf = 0.5 * facearea * (left.flux + right.flux - diss - (left.cons + right.cons) * vfn);
        auto nf = 0.5 * facearea * (left.flux + right.flux - diss);
        //for (int ii=0; ii<5; ++ii)
        //{
            //if (ii == 1) continue;
            //if (std::abs(nf[ii]) > 1e-5)
            //{
                //std::cout << ii << " " << left.flux[ii] << " " << right.flux[ii] << " " << diss[ii] << " " << facearea << " " << nf[ii] << " " << vfn << std::endl;
            //}
            //assert(std::abs(nf[ii]) < 1e-5);
        //}
        return nf;
    }

    auto RiemannSolver::left_eigenv(double gamma)
    {
        // inverse(R) or inverse(right_eigenv)

        //double nx = n(0);
        //double ny = n(1);
        //double nz = n(2);

        //double gamStar = gamma - 1.;

        //double gsk = 0.5 * gamStar * k / asq;
        //double gsu = 0.5 * gamStar * u / asq;
        //double gsv = 0.5 * gamStar * v / asq;
        //double gsw = 0.5 * gamStar * w / asq;
        //double gs = 0.5 * gamStar / asq;

        //double lstqn = 0.5 * a * qn / asq;
        //double lstnx = 0.5 * a * nx / asq;
        //double lstny = 0.5 * a * ny / asq;
        //double lstnz = 0.5 * a * nz / asq;

        Matrix5 L;

        double gc = gamma - 1.;

        auto vnet = u;

        L(0,0) = H + a * (vnet - a) / gc;
        L(0,1) = -(vnet + a / gc);
        L(0,2) = -v;
        L(0,3) = -w;
        L(0,4) = 1.;

        L(1,0) = -2. * H + 4 * a * a / gc;
        L(1,1) = 2. * vnet;
        L(1,2) = 2. * v;
        L(1,3) = 2. * w;
        L(1,4) = -2.;

        L(2,0) = -2. * v * a * a / gc;
        L(2,1) = 0.;
        L(2,2) = 2. * a * a / gc;
        L(2,3) = 0.;
        L(2,4) = 0.;

        L(3,0) = -2. * w * a * a / gc;
        L(3,1) = 0.;
        L(3,2) = 0.;
        L(3,3) = 2. * a * a / gc;
        L(3,4) = 0.;

        L(4,0) = H - a * (vnet + a) / gc;
        L(4,1) = -vnet + a / gc;
        L(4,2) = -v;
        L(4,3) = -w;
        L(4,4) = 1.;

        //

        /*L(0,0) = gsk + lstqn;
        L(0,1) = -gsu - lstnx;
        L(0,2) = -gsv - lstny;
        L(0,3) = -gsw - lstnz;
        L(0,4) = gs;

        L(1,0) = asq - 2. * gsk;
        L(1,1) = 2. * gsu;
        L(1,2) = 2. * gsv;
        L(1,3) = 2. * gsw;
        L(1,4) = -2. * gs;

        L(2,0) = gsk - lstqn;
        L(2,1) = -gsu + lstnx;
        L(2,2) = -gsv + lstny;
        L(2,3) = -gsw + lstnz;
        L(2,4) = gs;

        if (std::abs(nx) >= std::abs(ny) && std::abs(nx) >= std::abs(nz))
        {
            L(3,0) = (v - qn * ny) / nx;
            L(3,1) = ny;
            L(3,2) = (std::pow(ny,2) - 1.) / nx;
            L(3,3) = (ny * nz) / nx;
            L(3,4) = 0.;

            L(4,0) = (qn * nz - w) / nx;
            L(4,1) = -nz;
            L(4,2) = -(ny * nz) / nx;
            L(4,3) = (1. - std::pow(nz,2)) / nx;
            L(4,4) = 0.;
        }
        else if (std::abs(ny) >= std::abs(nx) && std::abs(ny) >= std::abs(nz))
        {
            L(3,0) = (qn * nx - u) / ny;
            L(3,1) = (1. - std::pow(nx,2)) / ny;
            L(3,2) = -nx;
            L(3,3) = -(nx * nz) / ny;
            L(3,4) = 0.;

            L(4,0) = (w - qn * nz) / ny;
            L(4,1) = (nx * nz) / ny;
            L(4,2) = nz;
            L(4,3) = (std::pow(nz,2) - 1.) / ny;
            L(4,4) = 0.;
        }
        else if (std::abs(nz) >= std::abs(nx) && std::abs(nz) >= std::abs(ny))
        {
            L(3,0) = (u - qn * nx) / nz;
            L(3,1) = (std::pow(nx,2) - 1.) / nz;
            L(3,2) = (nx * ny) / nz;
            L(3,3) = nx;
            L(3,4) = 0.;

            L(4,0) = (qn * ny - v) / nz;
            L(4,1) = -(nx * ny) / nz;
            L(4,2) = (1. - std::pow(ny,2)) / nz;
            L(4,3) = -ny;
            L(4,4) = 0.;
        }
        */
    
        return L;
    }

    Matrix5 RiemannSolver::Jacobian(const Matrix5& ws, const Matrix5& R, double gamma)
    {
        return ((R * ws) * left_eigenv(gamma));
    }

    void RiemannSolver::roe(const State& left, const State& right, Vector5& numflux, Matrix5& Aroe, double& max_eigen, double signed_area, double gamma)
    {
        calc_roe_ave_vars(left, right, gamma);

        Matrix5 R = right_eigenv();
        Matrix5 ws = abs_eigen(left, right);
        numflux = numerical_flux(left, right, ws, R, signed_area, gamma);
        Aroe = Jacobian(ws, R, gamma);
        max_eigen = max(ws);
        assert(max_eigen > 0.);

        //calc_roe_ave_vars(leftorig, rightorig, gamma, false);
        //R = right_eigenv(vfn);
        //ws = abs_eigen(vfn, leftorig, rightorig);
        //Aroe = Jacobian(ws, R, gamma, vfn);

        //if (abs(numflux(3)) >= 1e-0)
        //{
            //std::cout << numflux << std::endl;
            //std::cout << "LEFT VALUES" << std::endl;
            //left.print();
            //std::cout << "RIGHT VALUES" << std::endl;
            //right.print();
        //}
        
        //assert(abs(numflux(3)) < 1e-0);
    }

    std::tuple<double, double> RiemannSolver::speed_estimate_hllc(const State& left, const State& right, double gamma)
    {
        if (speed_estimate_hllc_ == SpeedEstimateHLLC::roe)
        {
            calc_roe_ave_vars(left, right, gamma); // needed if Roe vars are used to compute wave speed.
            return rhll_ws_est(left, right);
        }
        else if (speed_estimate_hllc_ == SpeedEstimateHLLC::pressure)
        {
            return rhll_ws_est_pres(left, right, gamma);
        }
    }

    std::tuple<Vector5, double, double, double> RiemannSolver::hllc(const State& left, const State& right, double signed_area, double gamma)
    {
        auto [sL, sR] = speed_estimate_hllc(left, right, gamma);

        Matrix5 ws = abs_eigen(left, right);
        auto max_eigen = max(ws);
        assert(max_eigen > 0.);
        auto numflux = hllc_numerical_flux(sL, sR, left, right, signed_area);

        return std::tie(numflux, max_eigen, sL, sR);

    }
}
