#pragma once
#include "../Force.hpp"
#include "../../data/Primitives.hpp"
#include "HeuresticBA.hpp"
#include "NativeBA.hpp"

namespace mdk {
    /**
     * A force field for bond angles (both heurestic and native variants,
     * depending on whether a triple has an associated native bond angle (in
     * such cases the native variant supercedes the heurestic variant).
     */
    class BondAngles: public Force {
    private:
        friend class HeuresticBA;
        /// The "heurestic" part of the potential.
        HeuresticBA* heurBA = nullptr;

        friend class NativeBA;
        /// The "native" part of the potential.
        NativeBA* natBA = nullptr;

        /// Whether a triple (i-1, i, i+1) is connected, i.e. in one chain.
        Bytes inRange;

    public:
        /**
         * Bind the force field to a simulation.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation &simulation) override;

        /**
         * "Asynchronous" part of the potential. Here it's the only part.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        void asyncPart(Dynamics &dynamics) override;
    };
}
