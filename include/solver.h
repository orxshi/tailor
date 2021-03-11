#ifndef SOLVER_H
#define	SOLVER_H

#include "read_mesh.h"
#include "limiter.h"
#include "riemann.h"
#include "bc.h"
#include "partition.h"
#include "var_exchanger.h"
#include "nei_exchanger.h"
#include "donor_info2.h"
#include "donor_request_exchanger.h"
#include "cell_exchanger.h"
#include "disconnecter.h"
#include "arrival_connecter.h"
//#include "facerequest_exchanger.h"
#include "mgmres.hpp"

namespace Tailor
{   
    class Solver
    {           
        public:     

            enum class RiemannSolverType
            {
                roe,
                hllc, // currently only for explicit formulation.
            };
    
            Solver(boost::mpi::communicator* comm, const std::vector<std::string>& filename, Profiler* profiler, Partition* partition=nullptr); 
            Solver();

            void set_partition(Partition* partition);

            void set_comm(boost::mpi::communicator* comm);
            void set_profiler(Profiler* prof);
                    
            void solve();
            double dt() const;
            void read_mesh(const std::vector<std::string>& file_name);
            const Partition* partition() const;
            Partition* partition();
            void transfer_oga_cell_type(const ArrCon<DonorInfo2>& di);
            void set_oga_cell_type_all_field();
            void rotate(const Tag& mesh, double ang, int axis, const vec3<double>& pivot);
            void exchange();
            void reset_oga_status();
            void reconnectivity();
            void update_fs(const Freestream& fs);
            bool repartition();
            void save_solution();
            const boost::mpi::communicator* comm() const;
            void read_settings();

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & dual_ts_;
                ar & repart_ratio_;
                ar & initratio_;
                ar & print_imbalance_;
                ar & print_vtk_;
                ar & print_repart_info_;
                ar & rebalance_thres_;
                ar & can_rebalance_;
                ar & force_rebalance_;
                ar & save_solution_;
                ar & restore_solution_;
                ar & half_cfl_;
                ar & print_map_;
                ar & fs_;
                ar & pseudo3D_;
                ar & omega_;
                ar & show_inner_res_;
                ar & show_inner_norm_;
                ar & print_outer_norm_;
                ar & sorder_;
                ar & torder_;
                ar & maxtimestep_;
                ar & tol_;
                ar & nsweep_;
                ar & cfl_;
                ar & delta_cfl_;
                ar & cfl_increase_freq_;
                ar & ncfl_increase_;
                ar & progressive_cfl_;
                ar & dt_; 
                ar & steady_;   
                ar & verbose_;
                ar & printfreq_;
                ar & finaltime_;
                ar & tempo_discre_;
                ar & nsolve_;
                ar & bc_;
                ar & load_estim_type_;
                ar & make_load_balance_;
                ar & init_max_res_;
                ar & last_max_res_;
                ar & riemann_solver_type_;
            }

        private:    
                
            bool dual_ts_;
            RiemannSolverType riemann_solver_type_;
            int repart_ratio_;
            double initratio_;
            bool print_imbalance_;
            bool print_vtk_;
            bool print_repart_info_;
            double rebalance_thres_;
            bool can_rebalance_;
            bool force_rebalance_;
            Profiler* profiler_;
            bool save_solution_;
            bool restore_solution_;
            bool half_cfl_;
            bool print_map_;
            Freestream fs_;
            bool pseudo3D_;
            double omega_;
            Gradient gradient_;
            bool show_inner_res_;
            bool show_inner_norm_;
            bool print_outer_norm_;
            int sorder_;
            int torder_;
            int maxtimestep_;
            double tol_;
            double nsweep_;
            double cfl_;
            double delta_cfl_;
            int cfl_increase_freq_;
            int ncfl_increase_;
            bool progressive_cfl_;
            double dt_; 
            bool steady_;   
            bool verbose_;
            int printfreq_;
            double finaltime_;
            std::string tempo_discre_;
            int nsolve_;
            BoundaryCondition bc_;

            Partition* partition_;
            boost::mpi::communicator* comm_;
            LoadEstim load_estim_type_;
            bool make_load_balance_;
            VarExchanger* var_exc_;
            VarExchanger* donor_var_exc_;

            double init_max_res_;
            double last_max_res_;
                
            void RK4(Mesh& mesh);
            void restore_solution();
            void connect_partition_cells();
            void exchange_ghosts();
            void solve_(int rank);
            void set_from_settings();
            //void calc_steady(Mesh& mesh, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts);
            void calc_steady(Mesh& mesh, int rank);
            void calc_steady();
            void update_cons_explicitly(Mesh& mesh);
            void update_cons_implicitly(Mesh& mesh, int ntimestep);
            void calc_R(Mesh& mesh);
            //void update_matrices(Mesh& mesh, const vararray& flux, const Matrix<NVAR, NVAR>& Aroe, MeshFace& mf, MeshCell& LC, MeshCell& RC, const vec3<double>& n, double vgn, double facearea);
            void update_matrices(Mesh& mesh, const vararray& flux, const Matrix<NVAR, NVAR>& Aroe, MeshFace* myface, MeshFace* commonface, MeshCell& LC, MeshCell& RC, double facearea, double volume, double vfn, const vec3<double>& vf, const State& left, const State& right, const varmat& TT);
            //void rhhl_update_matrices(const Matrix<NVAR, NVAR>& Aroe, MeshFace* myface, MeshFace* commonface, MeshCell& LC, MeshCell& RC, const vec3<double>& n, double vgn, double facearea, double SLm, double SRp, double vfn);
            //void hhl_update_matrices(MeshFace* myface, MeshFace* commonface, MeshCell& LC, MeshCell& RC, double facearea, double SLm, double SRp, const varmat& T, double vfn);
            void limit_prim_cons(Mesh& mesh, const MeshCell& mc, const MeshFace& mf, vararray& prim);
            bool sor(Mesh& mesh, int ntimestep);
            void tempo_discre(Mesh& mesh, bool calc_dtau);
            //void residual(Mesh& mesh);
            vararray resi(Mesh& mesh);
            void sweep(Mesh& mesh, MeshCell& mc, vararray& r, vararray& r2, vararray& r3, int& maxcell);
            Matrix<NVAR,NVAR> slipwall_M(const vec3<double>& n);
            void gmres(Mesh& mesh);
            void oga_interpolate(Mesh& mesh);
    };
}

#endif
