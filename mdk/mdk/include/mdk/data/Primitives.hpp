#pragma once
#include <Eigen/Dense>
#include "Vectors.hpp"

namespace mdk {
    /// A list of scalars (i.e. doubles).
    using Scalars = Eigen::VectorXd;

    /// A vector (writing \p Eigen::Vector3d everywhere is annoying).
    using Vector = Eigen::Vector3d;

    /**
     * Constant reference to a vector. Given the size of a Vector, passing
     * copies wouldn't be inefficient; it's here chiefly to hide the warnings.
     */
    using VRef = Vector const&;

    /**
     * List of small integers (usually we use this to store various enums
     * and true/false values).
     */
    using Bytes = std::vector<int8_t>;

    using Integers = std::vector<int>;

    /// List of pairs (used in storing VL pairs).
    using Pairs = std::vector<std::pair<int, int>>;
}
