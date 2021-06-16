#pragma once
#include "Primitives.hpp"
#include "../files/param/Parameters.hpp"
#include "../model/Model.hpp"

namespace mdk {
    /**
     * A data structure containing the charges of residues in a model, as
     * derived from a parameter file, in an array form.
     */
    class Charges: public Scalars {
    public:
        Charges() = default;
        Charges(Model const& model, param::Parameters const& params);
    };
}