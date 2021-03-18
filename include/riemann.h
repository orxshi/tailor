#ifndef POES_ROE_H
#define POES_ROE_H

#include "matrix.h"
#include "meshface.h"
#include "state.h"
#include "utility.h"

namespace Tailor
{
    enum class RiemannSolverType
    {
        roe,
        hllc
    };

    class RiemannSolver
    {
        public:

            RiemannSolver(RiemannSolverType riemann_solver_type, const State& left_state, const State& right_state, double face_area, double gamma, double& max_eigen, Vector5& flux);

            double u;
            double v;
            double w;
            double H;
            double k;
            double a;

        private:

            Vector5 hlle_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            Vector5 hllc_numerical_flux(double SLm, double SRp, const State& left, const State& right, double facearea);
            void rhll_ws_est(const State& left, const State& right, double& SLm, double& SRp);
            void rhll_ws_est_pres(const State& left, const State& right, double& SLm, double& SRp, double gamma);
            //void ws_est_pbased(const State& left, const State& right, double& SLm, double& SRp, double gamma);
            //Matrix5 rhll_abs_eigen(double alpha1, double alpha2, double SLm, double SRp, double vfn);
            void calc_roe_ave_vars(const State& left, const State& right, double gamma);
            Matrix5 right_eigenv();
            auto left_eigenv(double gamma);
            Matrix5 abs_eigen(const State& left, const State& right);
            Matrix5 eigen();
            Vector5 wave_strength_ver_1(const State& left, const State& right, double gamma);
            Vector5 wave_strength_ver_2(const State& left, const State& right, const Vector3& n);
            Vector5 dissipation_term(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double gamma);
            //Vector5 jump_ver_1(const State& left, const State& right);
            Vector5 jump_ver_2(const State& left, const State& right);
            //Vector5 jump_ver_3(const State& left, const State& right);
            Matrix5 Jacobian(const Matrix5& ws, const Matrix5& R, double gamma);
            Vector5 numerical_flux(const State& left, const State& right, const Matrix5& ws, const Matrix5& R, double facearea, double gamma);
            Vector5 rhll_numerical_flux(double SLm, double SRp, const State& left, const State& right, const State& leftn2, const State& rightn2, const Matrix5& ws, const Matrix5& R, double facearea, const Vector3& n, double gamma);
            void roe(const State& left, const State& right, Vector5& numflux, Matrix5& Aroe, double& max_eigen, double signed_area, double gamma);
            void hlle(const State& left, const State& right, Vector5& numflux, double& max_eigen, double signed_area, double gamma, double& SLm ,double& SRp);
            void hllc(const State& left, const State& right, Vector5& numflux, double& max_eigen, double signed_area, double gamma, double& SLm ,double& SRp);
    };
}

#endif
