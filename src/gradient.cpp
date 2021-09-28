#include "gradient.h"

namespace Tailor
{
    std::array<Vector3, NVAR> Gradient::ls_grad(const Mesh& mesh, const MeshCell& mc)
    {    
        std::array<Vector3, NVAR> data_;

        for (int i=0; i<NVAR; ++i)
        {
            data_[i] = 0.;
        }

        int ff=0;
        for (const MeshFace& f: mc.face())
        //for (const auto& f: mc.pnei())
        {
            //const auto& nei = mesh.cell(f);

            const MeshCell* nei = opposing_nei(mesh, f, mc.tag());
            assert(nei != nullptr);

            for (int i=0; i<NVAR; ++i)
            {
                //if (nei.btype() == BouType::partition) {
                    //data_[i] = 0.;
                    //continue;
                //}

                if (std::isnan(mc.ls_wx_[ff]))
                {
                    std::cout << "ls_wx size: " << mc.ls_wx_.size() << std::endl;
                    std::cout << "ff: " << ff << std::endl;
                    std::cout << "mc.pnei().size(): " << mc.pnei().size() << std::endl;
                }
                if (std::isnan(mc.ls_wy_[ff]))
                {
                    std::cout << "mesh: " << mesh.tag()() << std::endl;
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
                    std::cout << "ls_wy size: " << mc.ls_wy_.size() << std::endl;
                    std::cout << "ff: " << ff << std::endl;
                    std::cout << "mc.pnei().size(): " << mc.pnei().size() << std::endl;
                }
                assert(!std::isnan(mc.ls_wx_[ff]));
                assert(!std::isnan(mc.ls_wy_[ff]));
                assert(!std::isnan(mc.ls_wz_[ff]));
                //assert(!std::isnan(mesh.cell(f).prim(i)));
                assert(!std::isnan(mc.prim(i)));

                double tempf = nei->prim(i) - mc.prim(i);
                //double tempf = f.const_addr()->prim(i) - mc.prim(i);

                if (ff >= mc.ls_wx_.size())
                {
                    std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "ff: " << ff << std::endl;
                    std::cout << "mc.pnei().size(): " << mc.pnei().size() << std::endl;
                }

                assert(ff < mc.ls_wx_.size());
                assert(ff < mc.ls_wy_.size());
                assert(ff < mc.ls_wz_.size());

                data_[i](0) = data_[i](0) + mc.ls_wx_[ff] * tempf;
                data_[i](1) = data_[i](1) + mc.ls_wy_[ff] * tempf;
                data_[i](2) = data_[i](2) + mc.ls_wz_[ff] * tempf;
            }

            //if (nei.btype() == BouType::partition) {
                //break;
            //}

            ++ff;
        }

        for (int i=0; i<NVAR; ++i)
        {
            assert(!data_[i].isnan());
            //assert(!std::isnan(data_[i](0)));
            //assert(!std::isnan(data_[i](1)));
            //assert(!std::isnan(data_[i](2)));
        }

        return data_;
    }

    void Gradient::calc_ls_coef(Mesh& mesh)
    {
        for (MeshCell& mc: mesh.cell_)
        {
            mc.ls_wx_.reserve(mc.face().size());
            mc.ls_wy_.reserve(mc.face().size());
            mc.ls_wz_.reserve(mc.face().size());

            double dx, dy, dz, a1, a2, a3, psi;
            double r_11 = 0.;
            double r_12 = 0.;
            double r_13 = 0.;
            double r_22 = 0.;
            double r_23 = 0.;
            double r_33 = 0.;

            if (mc.oga_cell_type() == OGA_cell_type_t::non_resident || mc.oga_cell_type() == OGA_cell_type_t::ghost) {
                continue;
            }

            assert(mc.btype() == BouType::interior || mc.btype() == BouType::partition);

            if (mc.face().size() <= 1)
            {
                std::cout << "face size: " << mc.face().size() << std::endl;
            }
            //assert(mc.face().size() > 1);

            for (MeshFace& f: mc.face_p())
            {
                if (f.parent_cell().size() != 2)
                {
                    std::cout << "mesh: " << mesh.tag()() << std::endl;
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "mf: " << f.tag() << std::endl;
                    std::cout << "mc oga: " << static_cast<int>(mc.oga_cell_type()) << std::endl;
                    std::cout << "mc btype: " << static_cast<int>(mc.btype()) << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(f.btype()) << std::endl;
                    std::cout << "parent size: " << f.parent_cell().size()<< std::endl;
                }
                assert(f.parent_cell().size() == 2);
                const MeshCell* nei = opposing_nei(mesh, f, mc.tag());
                auto d = nei->poly().centroid() - mc.poly().centroid();

                dx = d(0);
                dy = d(1);
                dz = d(2);

                r_11 += dx * dx;
                r_12 += dx * dy;
                r_13 += dx * dz;
                r_22 += dy * dy;
                r_23 += dy * dz;
                r_33 += dz * dz;
            }

            r_11 = std::sqrt(r_11);
            r_12 /= r_11;
            r_13 /= r_11;
            r_22 = std::sqrt(r_22 - std::pow(r_12,2));
            r_23 = (r_23 - r_12 * r_13) / r_22;
            r_33 = std::sqrt(r_33 - std::pow(r_13,2.) - std::pow(r_23,2.));

            if (r_22 == 0.)
            {
                for (MeshFace& f: mc.face_p())
                {
                    const MeshCell* nei = opposing_nei(mesh, f, mc.tag());
                    auto d = nei->poly().centroid() - mc.poly().centroid();

                    dx = d(0);
                    dy = d(1);
                    dz = d(2);

                    std::cout << "dx: " << dx << std::endl;
                    std::cout << "dy: " << dy << std::endl;
                    std::cout << "dz: " << dz << std::endl;
                }
            }
            assert(r_22 != 0.);

            for (MeshFace& f: mc.face_p())
            //for (const auto& f: mc.pnei())
            {
                const MeshCell* nei = opposing_nei(mesh, f, mc.tag());
                auto d = nei->poly().centroid() - mc.poly().centroid();

                //auto d = mesh.cell(f).poly().centroid() - mc.poly().centroid();
                //auto d = f.const_addr()->poly().centroid() - mc.poly().centroid();

                dx = d(0);
                dy = d(1);
                dz = d(2);

                a1 = dx / std::pow(r_11,2.);
                a2 = (dy - dx * r_12 / r_11) / std::pow(r_22,2.);
                psi = (r_12 * r_23 - r_13 * r_22) / (r_11 * r_22);
                a3 = (dz - dy * r_23 / r_22 + psi * dx) / std::pow(r_33,2.);

                mc.ls_wx_.push_back(a1 - a2 * r_12 / r_11 + psi * a3);
                mc.ls_wy_.push_back(a2 - a3 * r_23 / r_22);
                mc.ls_wz_.push_back(a3);
            }

            for (const auto& a: mc.ls_wx_)
            {
                assert(!std::isnan(a));
            }
            for (const auto& a: mc.ls_wy_)
            {
                assert(!std::isnan(a));
            }
            for (const auto& a: mc.ls_wz_)
            {
                assert(!std::isnan(a));
            }
        }
    }
}
