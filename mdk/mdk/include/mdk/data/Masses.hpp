#pragma once
#include "Primitives.hpp"
#include "../model/Model.hpp"

namespace mdk {
    /**
     * Data class containing masses associated with residues, in an array form.
     */
    class Masses: public Scalars {
    public:
        Masses() = default;
        explicit Masses(Model const& model);
    };
}
