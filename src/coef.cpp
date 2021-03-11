#include "coef.h"

namespace Tailor
{
    /*void pres_coef(const Mesh& mesh, std::vector<double>& pc, double& lift_coef, const SpatialPartitionContainer* spc, const Freestream& fs)
    {
        // TODO need to insert Aref which is a function of Rref.
        const auto& wall = mesh.wall_boundaries();

        //double total_area = 0.;
        //for (const MeshCell& mc: wall) {
            //total_area += std::abs(mc.face()[0].face().signed_area());
        //}

        double pinf = fs.pinf_;
        double rhoinf = fs.rhoinf_;
        double gamma = fs.gamma_;
        double cinf = std::sqrt(gamma * pinf / rhoinf);

        //double mach_ref = fs.machair_ - fs.machfoil_;
        //double mach_ref = fs.machair_;
        //double dp; // dynamic pressure
        //double dp = 0.5 * fs.rhoinf_ * std::pow(mach_ref, 2.); // dynamic pressure
        double dp = 0.5 * fs.rhoinf_ * std::pow(96.443, 2.) * 6.406; // dynamic pressure
        //double dpa = dp * total_area;

        // total force on airfoil
        Vector3 F(0., 0., 0.);
        for (auto mc = wall.begin(); mc != wall.end(); ++mc)
        {
            const MeshCell& in = mesh.cell(mc->interior_boundary());

            assert(mc->face().size() == 1);
            const auto& mf = mc->face()[0];

            //double clocal = std::sqrt(gamma * mc.prim(4) / mc.prim(0));
            //double mach_ref = fs.machair_ - mf.vgn() / clocal;
            //double dp = 0.5 * fs.rhoinf_ * std::pow(mach_ref, 2.); // dynamic pressure

            Tag temp;
            if (!spc->is_resident(in.poly().centroid(), -1, temp)) {
                continue;
            }

            int i = std::distance(wall.begin(), mc);
            pc[i] = (mc->prim(4) - fs.pinf_) / dp;
            F = F + mc->face()[0].face().normal() * mc->prim(4) * std::abs(mc->face()[0].face().signed_area());
        }

        lift_coef = F(1) / dp;
    }

    void get_coef(const Mesh& mesh, int rank, const boost::mpi::communicator* comm, const SpatialPartitionContainer* spc, int iter, const Freestream& fs)
    {
        double local_lc;
        const auto& wall = mesh.wall_boundaries();
        std::vector<double> pc(wall.size());

        pres_coef(mesh, pc, local_lc, spc, fs);
        double global_lc = 0.;
        boost::mpi::all_reduce(*comm, local_lc, global_lc, std::plus<double>());

        {
            std::string fn = "pres_coef-";
            fn.append(std::to_string(rank));
            fn.append("-");
            fn.append(std::to_string(mesh.tag()()));
            fn.append("-");
            fn.append(std::to_string(iter));
            fn.append(".dat");

            std::ofstream out;
            out.open(fn);

            for (auto mc = wall.begin(); mc != wall.end(); ++mc)
            {
                out << mc->poly().centroid()(0);
                out << " "; 
                out << mc->poly().centroid()(1);
                out << " "; 
                int i = std::distance(wall.begin(), mc);
                out << -pc[i];
                out << std::endl;
            }

            out.close();
        }

        if (comm->rank() == 0)
        {
            std::string fn = "lift_coef-";
            fn.append(std::to_string(mesh.tag()()));
            fn.append("-");
            fn.append(std::to_string(iter));
            fn.append(".dat");

            std::ofstream out;
            out.open(fn, std::fstream::app);

            out << (fs.aoa_foil_x_ - 180.) + fs.aoa_air_x_;
            out << " "; 
            out << (90. - (fs.aoa_foil_x_ - 180.)) + fs.aoa_air_x_;
            out << " "; 
            out << (fs.aoa_foil_z_ - 180.) + fs.aoa_air_z_;
            out << " "; 
            out << global_lc; 
            out << std::endl; 
            out.close();
        }
    }

    void print_wall(const Mesh& mesh, int rank, int iter, const SpatialPartitionContainer* spc)
    {
        std::ofstream out;

        const auto& wall = mesh.wall_boundaries();
        if (wall.empty()) {
            return;
        }

        std::string fn = "wallcoef-";
        fn.append(std::to_string(rank));
        fn.append("-");
        fn.append(std::to_string(mesh.tag()()));
        fn.append("-");
        fn.append(std::to_string(iter));
        fn.append(".dat");
        out.open(fn);

        for (auto mc = wall.begin(); mc != wall.end(); ++mc)
        {
            const MeshCell& in = mesh.cell(mc->interior_boundary());
            const auto& mf = mc->face()[0];

            Tag temp;
            if (!spc->is_resident(in.poly().centroid(), -1, temp)) {
                continue;
            }

            out << mf.face().normal(0); out << " ";
            out << mf.face().normal(1); out << " ";
            out << mf.face().normal(2); out << " ";
            out << mf.face().signed_area(); out << " ";
            out << mf.face().centroid()(0); out << " ";
            out << mf.face().centroid()(1); out << " ";
            out << mf.face().centroid()(2); out << " ";
            out << mc->prim(4);
            out << std::endl;
        }

        out.close();

        {
            std::string fn = "wallpara-";
            fn.append(std::to_string(rank));
            fn.append("-");
            fn.append(std::to_string(mesh.tag()()));
            fn.append("-");
            fn.append(std::to_string(iter));
            fn.append(".vtk");
            mesh.print_wall_as_vtk(fn);
        }
    }*/
}
