#pragma once
#include "ESBase.hpp"

namespace mdk {
    /**
     * Debye-Hueckel screened electrostatic potential with a distance-dependent
     * electric permittivity (specifically eps_r = r_0/r).
     */
    class RelativeDH: public ESBase {
    protected:
        /**
         * Generate a VL spec.
         * @return Generated VL spec.
         */
        vl::Spec spec() const override;

    public:
        double screeningDist = 10.0 * angstrom;
        double r0 = 4.0 * angstrom;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        void asyncPart(Dynamics &dynamics) override;
    };
}
