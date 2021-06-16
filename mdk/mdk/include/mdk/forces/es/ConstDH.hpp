#pragma once
#include "ESBase.hpp"

namespace mdk {
    /**
     * Debye-Hueckel screeened electrostatic potential with constant
     * electric permittivity of the medium.
     */
    class ConstDH: public ESBase {
    protected:
        /**
         * Generate a VL spec.
         * @return Generated VL spec.
         */
        vl::Spec spec() const override;

    public:
        double screeningDist = 10.0 * angstrom;
        double permittivity = 80.0 * epsilon_0;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        void asyncPart(Dynamics &dynamics) override;
    };
}
