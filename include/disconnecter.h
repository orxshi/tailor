#ifndef DISCONNECTER_H
#define DISCONNECTER_H

#include "cell_exchanger.h"

namespace Tailor
{
    class Disconnecter
    {
        public:

            Disconnecter(SpatialPartitionContainer*, const Exchanger<Cell>*, Profiler*);
            void disconnect(bool mergebins, int rank);

        private:

            const Exchanger<Cell>* exc_;
            int cellsize(int fn) const;
            Profiler* profiler_;
            CellExchanger* cell_exchanger_;
            SpatialPartitionContainer* spc_;

            void disconnect_nonoverlapping_cells(bool mergebins, int rank);
            void disconnect_neiless_cells(int rank);
            void disconnect_empty_meshes();
            void disconnect_isolated_points();
            void disconnect_orphan_faces();
            void remove_ghosts();
    };
}

#endif
