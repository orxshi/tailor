#ifndef POES_ROE_H
#define POES_ROE_H

#include "matrix.h"
#include "meshface.h"
#include "state.h"
#include "utility.h"

namespace Tailor
{
    class RiemannSolver
    {
        public:

            RiemannSolver(const State& left, const State& right, double gamma, bool isbou);
            void roe(const State& left, const State& right, const State& leftorig, const State& rightorig, Vector5& numflux, Matrix5& Aroe, double& max_eigen, double signed_area, double gamma, double vfn);
            //void rhll(const State& left, const State& right, Vector5& numflux, Matrix5& Aroe, double& max_eigen, double signed_area, const Vector3& normal, double gamma, bool isbou, double& SLm ,double& SRp);
            void hlle(const State& left, const State& right, Vector5& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm ,double& SRp, double vfn);
            void hllc(const State& left, const State& right, Vector5& numflux, double& max_eigen, double signed_area, double gamma, bool isbou, double& SLm ,double& SRp, double vfn);

            double u;
            double v;
            double w;
            double H;
            double k;
            double a;

        private:

            Vector5 hlle_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            Vector5 hllc_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            void rhll_ws_est(const State& left, const State& right, double& SLm, double& SRp, double vfn);
            void rhll_ws_est_pres(const State& left, const State& right, double& SLm, double& SRp, double vfn, double gamma);
            //void ws_est_pbased(const State& left, const State& right, double& SLm, double& SRp, double gamma);
            //Matrix5 rhll_abs_eigen(double alpha1, double alpha2, double SLm, double SRp, double vfn);
            void calc_roe_ave_vars(const State& left, const State& right, double gamma, bool isbou);
            Matrix5 right_eigenv(double vfn);
            auto left_eigenv(double gamma, double vfn);
            Matrix5 abs_eigen(double vfn, const State& left, const State& right);
            Matrix5 eigen(double vfn);
            Vector5 wave_strength_ver_1(const State& left, const State& right, double gamma, double vfn);
            Vector5 wave_strength_ver_2(const State& left, const State& right, const Vector3& n);
            Vector5 dissipation_term(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double gamma, double vfn);
            //Vector5 jump_ver_1(const State& left, const State& right);
            Vector5 jump_ver_2(const State& left, const State& right);
            //Vector5 jump_ver_3(const State& left, const State& right);
            Matrix5 Jacobian(const Matrix5& ws, const Matrix5& R, double gamma, double vfn);
            Vector5 numerical_flux(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double facearea, double gamma, double vfn);
            Vector5 rhll_numerical_flux(double SLm, double SRp, const State& left, const State& right, const State& leftn2, const State& rightn2, const Matrix5& ws, const Matrix5& R, double facearea, const Vector3& n, double gamma, double vfn);

            //double u;
            //double v;
            //double w;
            //double H;
            //double k;
            //double a;
    };
}

#endif
