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

#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/adapter/zero_copy.hpp>

namespace Tailor
{   
    enum class TemporalDiscretization
    {
        undefined = -1,
        euler = 0,
        three_time_level = 1,
        rk4 = 2,
    };

    enum class SpatialDiscretization
    {
        undefined = -1,
        upwind = 0,
        MUSCL = 1,
    };

    typedef amgcl::backend::builtin<double> Backend;

    //typedef amgcl::make_solver<
    //    amgcl::amg<
    //    Backend,
    //    amgcl::runtime::coarsening::wrapper,
    //    amgcl::runtime::relaxation::wrapper
    //        >,
    //    amgcl::runtime::solver::wrapper<Backend>
    //        > AMGCLSolver;



    //typedef amgcl::make_solver<
    //amgcl::preconditioner::dummy<Backend>
    //,
    //amgcl::solver::gmres<Backend>
    //>
    //AMGCLSolver;

    class Solver
    {           
        public:     
    
            Solver(boost::mpi::communicator* comm, const std::vector<std::string>& filename, Profiler* profiler, Partition* partition=nullptr); 
            Solver();

            void update_flow_field_after_mesh_motion();

            int nsolve() const;
            void print_settings() const;
            void set_partition(Partition* partition);

            void set_comm(boost::mpi::communicator* comm);
            void set_profiler(Profiler* prof);
                    
            void solve(int max_time_step);
            double dt() const;
            void read_mesh(const std::vector<std::string>& file_name);
            const Partition* partition() const;
            Partition* partition();
            void transfer_oga_cell_type(const ArrCon<DonorInfo2>& di);
            void set_oga_cell_type_all_field();
            void rotate(const Tag& mesh, double ang, const Vector3& axis, const Vector3& pivot);
            void move(const Tag& mesh, const Vector3& v);
            void exchange();
            void reset_oga_status();
            void reconnectivity();
            void update_fs(const Freestream& fs);
            Freestream fs() const;
            bool repartition();
            const boost::mpi::communicator* comm() const;
            void read_settings();

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & implicit_;
                ar & partition_;
                ar & limiter_type_;
                ar & dual_ts_;
                ar & repart_ratio_;
                ar & initratio_;
                ar & print_imbalance_;
                ar & print_vtk_init_;
                ar & print_repart_info_;
                ar & rebalance_thres_;
                ar & can_rebalance_;
                ar & force_rebalance_;
                ar & half_cfl_;
                ar & print_map_;
                ar & fs_;
                ar & pseudo3D_;
                ar & omega_;
                ar & show_inner_res_;
                ar & show_inner_norm_;
                ar & print_residual_;
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
                ar & temporal_discretization_;
                ar & nsolve_;
                ar & bc_;
                ar & load_estim_type_;
                ar & make_load_balance_;
                ar & initial_global_residual_;
                ar & last_global_residual_;
                ar & riemann_solver_type_;
                ar & increase_cfl_;
                ar & cfl_multiplier_;
                ar & cfl_ratio_;
                ar & global_nmesh_;
                ar & print_linear_solver_error_;
                ar & linear_solver_max_iteration_;
                ar & linear_solver_max_restart_;
                ar & linear_solver_abs_error_;
                ar & linear_solver_rel_error_;
                ar & print_vtk_only_last_step_;
                ar & print_vtk_interval_;
                ar & temporal_discretization_;
                ar & spatial_discretization_;
            }

        private:    

            bool implicit_;
            int print_vtk_interval_;
            LimiterType limiter_type_;
            bool print_vtk_only_last_step_;
            bool print_linear_solver_error_;
            int linear_solver_max_iteration_;
            int linear_solver_max_restart_;
            double linear_solver_abs_error_;
            double linear_solver_rel_error_;
            int global_nmesh_;
            double cfl_multiplier_;
            double cfl_ratio_;
            bool increase_cfl_;
            bool dual_ts_;
            RiemannSolverType riemann_solver_type_;
            int repart_ratio_;
            double initratio_;
            bool print_imbalance_;
            bool print_vtk_init_;
            bool print_repart_info_;
            double rebalance_thres_;
            bool can_rebalance_;
            bool force_rebalance_;
            Profiler* profiler_;
            //bool save_solution_;
            //bool restore_solution_;
            bool half_cfl_;
            bool print_map_;
            Freestream fs_;
            bool pseudo3D_;
            double omega_;
            bool show_inner_res_;
            bool show_inner_norm_;
            bool print_residual_;
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
            TemporalDiscretization temporal_discretization_;
            SpatialDiscretization spatial_discretization_;
            int nsolve_;
            BoundaryCondition bc_;
            Vector5 initial_global_residual_;
            Vector5 last_global_residual_;
            bool use_local_time_step_;

            Partition* partition_;
            boost::mpi::communicator* comm_;
            LoadEstim load_estim_type_;
            bool make_load_balance_;
            VarExchanger* var_exc_;
            VarExchanger* donor_var_exc_;
                
            Vector5 non_linear_iteration();
            void init_old_conservative_var();
            void compute_gradient_coef();
            void calc_mesh_velocities();
            void reset_overset_mesh_exchanger();
            void reset_partitioned_mesh_exchanger();
            void init_partitioned_mesh_exchanger();
            void linear_solver(Mesh &mesh);
            //void RK4(Mesh& mesh);
            void connect_partition_cells();
            void exchange_ghosts();
            void set_from_settings();
            //void calc_steady(Mesh& mesh, std::function<void(CellExchanger& cell_exchanger)> exchange_ghosts);
            //void update_cons_explicitly(Mesh& mesh);
            //void update_cons_implicitly(Mesh& mesh, int ntimestep);
            void calc_change_in_conserved_var(Mesh& mesh, int runge_kutta_stage);
            void evolve_solution_in_time(Mesh& mesh);
            void evolve_old_solution_in_time(Mesh& mesh);
            void compute_sum_of_fluxes(Mesh &mesh);
            void compute_sum_of_fluxes(Mesh& mesh, int ntimestep);
            bool sor(Mesh& mesh, int ntimestep);
            void temporal_discretization(Mesh& mesh);
            Vector5 compute_residual(Mesh& mesh);
            void sweep(Mesh& mesh, MeshCell& mc, Vector5& r, Vector5& r2, Vector5& r3, int& maxcell);
            Matrix5 slipwall_M(const Vector3& n);
            void gmres(Mesh& mesh);
            void oga_interpolate(Mesh& mesh);
            void update_matrices(MeshFace *this_face, MeshFace *common_face, MeshCell& left_cell, MeshCell& right_cell, double facearea, const Vector3& face_velocity, double gamma, const Matrix5& rotation_matrix, const Matrix5& inv_rotation_matrix, const Matrix5& Aroe, const State& left_state, const State& right_state);
            Vector5 face_conservative_var(const Mesh &mesh, const MeshCell &mc, const MeshFace &mf);
            void print_residual(const Vector5& residual);
            void print_mesh_vtk(std::string);
            void update_ghosts();
            void update_donors(Mesh& mesh);
            void set_boundary_conditions(Mesh& mesh);
            Vector5 get_global_residual(const Vector5& local_residual, int ntimestep);
            void increase_cfl(const Vector5& global_residual, int timestep);
            void print_sub_solver_residual(int ntimestep, const Vector5& residual);
            std::vector<double> amgcl(int n, int nz_num, const std::vector<int>& ia, const std::vector<int>& ja, const std::vector<double>& a, const std::vector<double>& rhs);
            std::vector<double> gmres(int n, int nz_num, std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a, std::vector<double>& rhs);
            void update_partitioned_mesh_exchanger();
            void update_overset_mesh_exchanger();
            void first_order_residual(Vector5& res, const MeshCell& mc);
            void compute_gradient(Mesh& mesh);
            void rotate_velocity(const Tag &mesh, double ang, int axis, const Vector3 &pivot);
            void print_status(int linear_iter, const Vector5& global_residual);
    };

    std::tuple<Matrix5, Matrix5> get_rotation_matrix(const Vector3& normal);
    std::tuple<Vector5, double, Matrix5> compute_flux(RiemannSolverType riemann_solver_type, double face_area, const Matrix5& inverse_rotation_matrix, const State& rotated_left_state, const State& rotated_right_state, const State& left_state, const State& right_state, double gamma);

    void parallel_print(const boost::mpi::communicator& comm, std::string message);
}

#endif
