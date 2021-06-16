#pragma once
#include "../Force.hpp"
#include <Eigen/Geometry>

namespace mdk {
    /**
     * Solid walls.
     */
    class SolidWall: public Force {
    public:
        Eigen::Hyperplane<double, 3> wall;

        void asyncPart(Dynamics& dynamics) override;
    };
}