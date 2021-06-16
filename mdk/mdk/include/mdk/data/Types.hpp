#pragma once
#include "Primitives.hpp"
#include "../model/Model.hpp"

namespace mdk {
    /**
     * Data class containing types of residues, in an array form.
     */
    class Types: public Eigen::Matrix<ResType, Eigen::Dynamic, 1> {
    public:
        Types() = default;
        explicit Types(Model const& model);
    };
}
