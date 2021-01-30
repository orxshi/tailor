#include "roe.h"

namespace Tailor
{
    Roe::Roe(const State& left, const State& right, double gamma, bool isbou)
    {
        calc_roe_ave_vars(left, right, gamma, isbou);
    }

    //vararray Roe::jump_ver_1(const State& left, const State& right)
    //{
    //    double rhoL = left.rho;
    //    double pL = left.p;
    //    //double qnL = left.qn;
    //    //double qlL = left.ql;
    //    //double qmL = left.qm;
    //    double uL = left.u;
    //    double vL = left.v;
    //    double wL = left.w;

    //    double rhoR = right.rho;
    //    double pR = right.p;
    //    //double qnR = right.qn;
    //    //double qlR = right.ql;
    //    //double qmR = right.qm;
    //    double uR = right.u;
    //    double vR = right.v;
    //    double wR = right.w;

    //    assert(!std::isnan(left.ql));
    //    assert(!std::isnan(right.ql));

    //    vararray jump;

    //    jump[0] = rhoR - rhoL;
    //    jump[1] = pR - pL;
    //    //jump[2] = qnR - qnL;
    //    //jump[3] = qlR - qlL;
    //    //jump[4] = qmR - qmL;
    //    jump[2] = uR - uL;
    //    jump[3] = vR - vL;
    //    jump[4] = wR - wL;

    //    return jump;
    //}
    
    vararray Roe::jump_ver_2(const State& left, const State& right)
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

        vararray jump;

        jump[0] = rhoR - rhoL;
        jump[1] = rhoR*uR - rhoL*uL; 
        jump[2] = rhoR*vR - rhoL*vL; 
        jump[3] = rhoR*wR - rhoL*wL; 
        jump[4] = ER - EL; 

        return jump;
    }

    //vararray Roe::jump_ver_3(const State& left, const State& right)
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

    //    vararray jump;

    //    jump[0] = rhoR - rhoL;
    //    jump[1] = uR - uL;
    //    jump[2] = vR - vL;
    //    jump[3] = wR - wL;
    //    jump[4] = pR - pL;

    //    return jump;
    //}

    void Roe::calc_roe_ave_vars(const State& left, const State& right, double gamma, bool isbou)
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
            std::cout << "isbou: " << isbou << std::endl;
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

    varmat Roe::right_eigenv(double vfn)
    {
        varmat R;

        double vnet = u - vfn;

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

    void rhll_normals(const State& left, const State& right, const vec3<double>& normal, double& alpha1, double& alpha2, vec3<double>& n1, vec3<double>& n2, bool isbou)
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

        n1.set(nx1, ny1, nz1);

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
            std::cout << "isbou: " << isbou << std::endl;
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

        n2.set(nx2, ny2, nz2);

        assert(!isnan(n1(0)));
        assert(!isnan(n1(1)));
        assert(!isnan(n1(2)));
        assert(!isnan(n2(0)));
        assert(!isnan(n2(1)));
        assert(!isnan(n2(2)));
    }

    //void Roe::ws_est_pbased(const State& left, const State& right, double& SLm, double& SRp, double gamma)
    //{
    //    double aL = left.a;
    //    double pL = left.p;
    //    double rhoL = left.rho;
    //    double qnL = left.qn;

    //    double aR = right.a;
    //    double pR = right.p;
    //    double rhoR = right.rho;
    //    double qnR = right.qn;

    //    auto ppvrs = 0.5 * (pL + pR) - 0.5 * (qnR - qnL) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
    //    auto ps = std::max(0., ppvrs);

    //    double q;
    //    if (ps <= pL)
    //    {
    //        q = 1.;
    //    }
    //    else
    //    {
    //        q = std::sqrt(1. + (gamma + 1) * (ps/pL - 1) / (2. * gamma));
    //    }

    //    SLm = qnL - aL * q;

    //    if (ps <= pR)
    //    {
    //        q = 1.;
    //    }
    //    else
    //    {
    //        q = std::sqrt(1. + (gamma + 1) * (ps/pR - 1) / (2. * gamma));
    //    }

    //    SRp = qnR + aR * q;
    //}

    void Roe::rhll_ws_est(const State& left, const State& right, double& SLm, double& SRp, double vfn)
    {
        double uL = left.u;
        double vL = left.v;
        double wL = left.w;
        double aL = left.a;

        double uR = right.u;
        double vR = right.v;
        double wR = right.w;
        double aR = right.a;

        SLm = std::min({0., u - a + vfn, uL - aL + vfn});
        SRp = std::max({0., u + a + vfn, uR + aR + vfn});

        SLm = std::min({0., u - a, uL - aL});
        SRp = std::max({0., u + a, uR + aR});
    }

    /*varmat Roe::rhll_abs_eigen(double alpha1, double alpha2, double SLm, double SRp, double vfn)
    {
        auto ws = abs_eigen(vfn);
        auto eig = eigen(vfn);
        varmat I;
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

    varmat Roe::eigen(double vfn)
    {
        // wave speeds

        varmat ws;
        ws = 0.;

        assert(!std::isnan(a));

        double vnet = u - vfn;

        ws(0,0) = vnet - a;
        ws(1,1) = vnet;
        ws(2,2) = vnet;
        ws(3,3) = vnet;
        ws(4,4) = vnet + a;

        return ws;
    }

    varmat Roe::abs_eigen(double vfn, const State& left, const State& right)
    {
        // Absolute values of the wave speeds

        varmat ws;
        ws = 0.;

        assert(!std::isnan(a));

        double vnet = u - vfn;

        ws(0,0) = std::abs(vnet - a);
        ws(1,1) = std::abs(vnet);
        ws(2,2) = std::abs(vnet);
        ws(3,3) = std::abs(vnet);
        ws(4,4) = std::abs(vnet + a);

        double lL[5];
        double lR[5];
        double l[5];

        lL[0] = left.u - left.a;
        lL[1] = left.u;
        lL[2] = left.u;
        lL[3] = left.u;
        lL[4] = left.u + left.a;

        lR[0] = right.u - right.a;
        lR[1] = right.u;
        lR[2] = right.u;
        lR[3] = right.u;
        lR[4] = right.u + right.a;

        l[0] = u - a;
        l[1] = u;
        l[2] = u;
        l[3] = u;
        l[4] = u + a;

        for (int i=0; i<5; ++i)
        {
            double sigma = std::max(0., l[i] - lL[i]);
            sigma = std::max(sigma, lR[i] - l[i]);

            double akL = l[i] - sigma; 
            double akR = l[i] + sigma; 

            //if (lL[i] < 0. && lR[i] > 0.)
            if (akL < 0. && akR > 0.)
            {
                //ws(i,i) = ((lR[i] + lL[i]) * l[i] - 2. * lR[i] * lL[i]) / (lR[i] - lL[i]);
                ws(i,i) = ((akR + akL) * l[i] - 2. * akR * akL) / (akR - akL);
            }
        }

        return ws;
    }

    /*vararray Roe::wave_strength_ver_2(const State& left, const State& right, const vec3<double>& n)
    {
        auto jump = jump_ver_2(left, right);

        vararray ldu;

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

    vararray Roe::wave_strength_ver_1(const State& left, const State& right, double gamma, double vfn)
    {
        vararray jump = jump_ver_2(left, right);
        vararray ldu;
        
        auto vnet = u - vfn;

        ldu[2] = jump[2] - v * jump[0];
        ldu[3] = jump[3] - w * jump[0];

        double j4 = jump[4] - ldu[2] * v - ldu[3] * w;

        ldu[1] = (gamma - 1.) * (jump[0] * (H - std::pow(vnet, 2.)) + vnet * jump[1] - j4) / std::pow(a, 2.); 
        ldu[0] = (jump[0] * (vnet + a) - jump[1] - a * ldu[1]) / (2. * a); 
        ldu[4] = jump[0] - (ldu[0] + ldu[1]);

        assert(!std::isnan(jump[0]));
        assert(!std::isnan(jump[1]));
        assert(!std::isnan(jump[2]));
        assert(!std::isnan(jump[3]));
        assert(!std::isnan(jump[4]));
        

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

    vararray Roe::dissipation_term(const State& left, const State& right, const varmat& ws, const varmat& R, double gamma, double vfn)
    {
        vararray diss;

        vararray ldu = wave_strength_ver_1(left, right, gamma, vfn);
        //auto ldu = wave_strength_ver_2(left, right, n);


        for (int i=0; i<NVAR; ++i)
        {
            diss[i] = 0.;
            for (int j=0; j<NVAR; ++j)
            {
                assert(!std::isnan(ws(j,j)));
                assert(!std::isnan(ldu[j]));
                assert(!std::isnan(R(i,j)));
                diss[i] += ws(j,j) * ldu[j] * R(i,j);
            }
        }

        return diss;
    }

    vararray Roe::rhll_numerical_flux(double SLm, double SRp, const State& left, const State& right, const State& leftn2, const State& rightn2, const varmat& ws, const varmat& R, double facearea, const vec3<double>& n, double gamma, double vfn)
    {
        vararray diss = dissipation_term(leftn2, rightn2, ws, R, gamma, vfn);
        auto nf = (SRp*left.flux - SLm*right.flux)/(SRp-SLm) - 0.5 * diss;
        nf = facearea * nf;
        return nf;
    }

    vararray Roe::hlle_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea)
    {
        auto nf = (SRp*left.flux - SLm*right.flux + SLm*SRp*(right.cons-left.cons))/(SRp-SLm);
        nf = facearea * nf;
        return nf;
    }

    vararray Roe::hllc_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea)
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

        vararray vecL;
        vecL[0] = 1.;
        vecL[1] = Ss;
        vecL[2] = vL;
        vecL[3] = wL;
        vecL[4] = EL / rhoL + (Ss - uL) * (Ss + pL / (rhoL * (SLm - uL)));

        vararray vecR;
        vecR[0] = 1.;
        vecR[1] = Ss;
        vecR[2] = vR;
        vecR[3] = wR;
        vecR[4] = ER / rhoR + (Ss - uR) * (Ss + pR / (rhoR * (SRp - uR)));

        auto usL = vecL * rhoL * ((SLm - uL) / (SLm - Ss));
        auto usR = vecR * rhoR * ((SRp - uR) / (SRp - Ss));

        auto fsL = fL + (usL - left.cons) * SLm;
        auto fsR = fR + (usR - right.cons) * SRp;

        if (Ss >= 0.)
        {
            return facearea * fsL;
        }
        else
        {
            return facearea * fsR;
        }
    }

    vararray Roe::numerical_flux(const State& left, const State& right, const varmat& ws, const varmat& R, double facearea, double gamma, double vfn)
    {
        vararray diss = dissipation_term(left, right, ws, R, gamma, vfn);
        assert(!std::isnan(left.flux[0]));
        assert(!std::isnan(right.flux[0]));
        assert(!std::isnan(diss[0]));
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

    auto Roe::left_eigenv(double gamma, double vfn)
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

        varmat L;

        double gc = gamma - 1.;

        auto vnet = u - vfn;

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

    varmat Roe::Jacobian(const varmat& ws, const varmat& R, double gamma, double vfn)
    {
        return ((R * ws) * left_eigenv(gamma, vfn));
    }

    void Roe::bbb(const State& left, const State& right, const State& leftorig, const State& rightorig, vararray& numflux, varmat& Aroe, double& max_eigen, double signed_area, double gamma, double vfn)
    {
        varmat R = right_eigenv(vfn);
        varmat ws = abs_eigen(vfn, left, right);
        numflux = numerical_flux(left, right, ws, R, signed_area, gamma, vfn);
        Aroe = Jacobian(ws, R, gamma, vfn);
        max_eigen = ws.max();
        assert(max_eigen > 0.);

        calc_roe_ave_vars(leftorig, rightorig, gamma, false);
        R = right_eigenv(vfn);
        ws = abs_eigen(vfn, leftorig, rightorig);
        Aroe = Jacobian(ws, R, gamma, vfn);

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

    //void Roe::rhll(const State& left, const State& right, vararray& numflux, varmat& Aroe, double& max_eigen, double signed_area, const vec3<double>& normal, double gamma, bool isbou, double& SLm ,double& SRp)
    //{
    //    if (isbou)
    //    {
    //        return bbb(left, right, numflux, Aroe, max_eigen, signed_area, normal, gamma);
    //    }
    //    double alpha1, alpha2;
    //    vec3<double> n1, n2;

    //    rhll_normals(left, right, normal, alpha1, alpha2, n1, n2, isbou);
    //    rhll_ws_est(left, right, SLm, SRp);

    //    State leftn2 = left;
    //    State rightn2 = right;

    //    leftn2.nx = n2(0); 
    //    leftn2.ny = n2(1); 
    //    leftn2.nz = n2(2); 

    //    rightn2.nx = n2(0); 
    //    rightn2.ny = n2(1); 
    //    rightn2.nz = n2(2); 

    //    calc_roe_ave_vars(leftn2, rightn2, gamma, isbou);
    //    varmat R = right_eigenv();
    //    auto ws = rhll_abs_eigen(alpha1, alpha2, SLm, SRp);
    //    numflux = rhll_numerical_flux(SLm, SRp, left, right, leftn2, rightn2, ws, R, signed_area, n2, gamma);

    //    Aroe = Jacobian(ws, R, n2, gamma);
    //    max_eigen = std::max({std::abs(ws(0,0)), std::abs(ws(1,1)), std::abs(ws(2,2)), std::abs(ws(3,3)), std::abs(ws(4,4))});
    //}

    void Roe::hlle(const State& left, const State& right, vararray& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm ,double& SRp, double vfn)
    {
        rhll_ws_est(left, right, SLm, SRp, vfn);
        //ws_est_pbased(left, right, SLm, SRp, gamma);

        varmat ws = abs_eigen(vfn, left, right);
        numflux = hlle_numerical_flux(SLm, SRp, left, right, signed_area);

        max_eigen = std::max({std::abs(ws(0,0)), std::abs(ws(1,1)), std::abs(ws(2,2)), std::abs(ws(3,3)), std::abs(ws(4,4))});
    }

    void Roe::hllc(const State& left, const State& right, vararray& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm , double& SRp, double vfn)
    {
        rhll_ws_est(left, right, SLm, SRp, vfn);
        //ws_est_pbased(left, right, SLm, SRp, gamma);

        varmat ws = abs_eigen(vfn, left, right);
        numflux = hllc_numerical_flux(SLm, SRp, left, right, signed_area);

        max_eigen = std::max({std::abs(ws(0,0)), std::abs(ws(1,1)), std::abs(ws(2,2)), std::abs(ws(3,3)), std::abs(ws(4,4))});
    }
}
