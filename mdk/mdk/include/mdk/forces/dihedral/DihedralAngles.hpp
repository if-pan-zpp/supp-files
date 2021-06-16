#pragma once
#include "../Force.hpp"
#include "../../data/Primitives.hpp"
#include "ComplexNativeDihedral.hpp"
#include "SimpleNativeDihedral.hpp"
#include "HeuresticDihedral.hpp"
#include <variant>

namespace mdk {
    /**
     * A force field for dihedral angles, composed of:
     * - the "native" part, which exists in two mutually exclusive variants
     *   (simple and complex);
     * - the "heurestic" part.
     * The native part supercedes heurestic part whenever native dihedral angle
     * is defined for a quadruple.
     */
    class DihedralAngles: public Force {
    private:
        friend class ComplexNativeDihedral;
        friend class SimpleNativeDihedral;
        /**
         * The "native" part, i.e. either nothing, complex variant or the
         * simple variant.
         */
        std::variant<std::monostate, ComplexNativeDihedral*, SimpleNativeDihedral*>
            natDih = std::monostate();

        friend class HeuresticDihedral;
        /// The "heurestic" part of the potential.
        HeuresticDihedral *heurDih = nullptr;

        Bytes inRange;

    public:
        /**
         * Bind the force field to a simulation.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * "Asynchronous" part of the potential. Here it's the only part.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        void asyncPart(Dynamics &dynamics) override;
    };
}
