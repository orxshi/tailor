#ifndef POES_ROE_H
#define POES_ROE_H

#include "matrix.h"
#include "meshface.h"
#include "state.h"
#include "utility.h"

namespace Tailor
{
    class Roe
    {
        public:

            Roe(const State& left, const State& right, double gamma, bool isbou);
            void bbb(const State& left, const State& right, const State& leftorig, const State& rightorig, vararray& numflux, varmat& Aroe, double& max_eigen, double signed_area, double gamma, double vfn);
            //void rhll(const State& left, const State& right, vararray& numflux, varmat& Aroe, double& max_eigen, double signed_area, const vec3<double>& normal, double gamma, bool isbou, double& SLm ,double& SRp);
            void hlle(const State& left, const State& right, vararray& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm ,double& SRp, double vfn);
            void hllc(const State& left, const State& right, vararray& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm ,double& SRp, double vfn);

            double u;
            double v;
            double w;
            double H;
            double k;
            double a;

        private:

            vararray hlle_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            vararray hllc_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            void rhll_ws_est(const State& left, const State& right, double& SLm, double& SRp, double vfn);
            //void ws_est_pbased(const State& left, const State& right, double& SLm, double& SRp, double gamma);
            //varmat rhll_abs_eigen(double alpha1, double alpha2, double SLm, double SRp, double vfn);
            void calc_roe_ave_vars(const State& left, const State& right, double gamma, bool isbou);
            varmat right_eigenv(double vfn);
            auto left_eigenv(double gamma, double vfn);
            varmat abs_eigen(double vfn, const State& left, const State& right);
            varmat eigen(double vfn);
            vararray wave_strength_ver_1(const State& left, const State& right, double gamma, double vfn);
            vararray wave_strength_ver_2(const State& left, const State& right, const vec3<double>& n);
            vararray dissipation_term(const State& left, const State& right, const varmat& ws, const varmat& R, double gamma, double vfn);
            //vararray jump_ver_1(const State& left, const State& right);
            vararray jump_ver_2(const State& left, const State& right);
            //vararray jump_ver_3(const State& left, const State& right);
            varmat Jacobian(const varmat& ws, const varmat& R, double gamma, double vfn);
            vararray numerical_flux(const State& left, const State& right, const varmat& ws, const varmat& R, double facearea, double gamma, double vfn);
            vararray rhll_numerical_flux(double SLm, double SRp, const State& left, const State& right, const State& leftn2, const State& rightn2, const varmat& ws, const varmat& R, double facearea, const vec3<double>& n, double gamma, double vfn);

            //double u;
            //double v;
            //double w;
            //double H;
            //double k;
            //double a;
    };
}

#endif
