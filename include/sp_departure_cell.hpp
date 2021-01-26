#ifndef SP_DEPARTURE_CELL_HPP
#define	SP_DEPARTURE_CELL_HPP

#include "tag.h"

namespace Tailor
{
    class DepartureCell
    {
        Tag icell_;
        Tag imesh_;

        public:

        DepartureCell() = default;
        DepartureCell(const Tag& icell, const Tag& imesh): icell_(icell), imesh_(imesh) {}

        const Tag& cell_tag() const;
        const Tag& mesh_tag() const;
        void set_cell_tag(const Tag& t);
        void set_mesh_tag(const Tag& t);
    };
}

#endif
