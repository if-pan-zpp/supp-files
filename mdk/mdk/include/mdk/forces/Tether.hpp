#pragma once
#include "Force.hpp"
#include "../kernels/Harmonic.hpp"
#include "../data/Chains.hpp"

namespace mdk {
    /**
     * Harmonic tether forces between consecutive residues in a chain.
     */
    class Tether: public Force {
    public:
        /**
         * Construct a \p Tether object.
         * @param fromNative Whether the equilibrium bond distances should be
         * derived from the native structure or set to a default value
         * (3.8 angstrem).
         */
        explicit Tether(bool fromNative);

        /**
         * Bind the object to the simulation. Additionally we fetch the
         * chain data to determine which of the pairs (i, i+1) are connected.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces
         * to.
         */
        void asyncPart(Dynamics &dynamics) override;

    private:
        /**
         * The underlying harmonic force kernel. It is separated from this class
         * because it's used in some other places, notably for standard SSBOND
         * interactions.
         */
        Harmonic harm;

        /**
         * Whether the equilibrium bond distances should be
         * derived from the native structure or set to a default value
         * (3.8 angstrem).
         */
        bool fromNative;

        /**
         * Number of residues, saved for convenience.
         */
        int n;

        /**
         * A list of equilibrium distances.
         */
        Scalars dist0;

        /**
         * isConnected[i] = 1 if pair (i, i+1) is entirely in a chain, i.e.
         * is tethered; 0 otherwise.
         */
        Bytes isConnected;
    };
}
