#pragma once
#include "NonlocalForce.hpp"
#include "../kernels/ShiftedTruncatedLJ.hpp"
#include "../data/Chains.hpp"

namespace mdk {
    /**
     * A "Pauli exclusion force" that keeps the residues from overlapping.
     */
    class PauliExclusion: public NonlocalForce {
    public:
        /**
         * A shifted and truncated version of the Lennard-Jones potential
         * that is used as the actual force.
         */
        ShiftedTruncatedLJ stlj;

        PauliExclusion();

        /**
         * Bind the object to the simulation.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and
         * forces to.
         */
        void asyncPart(Dynamics &dynamics) override;

        /**
         * An action performed when a Verlet list is reconstructed; here we
         * simply copy the list verbatim.
         */
        void vlUpdateHook() override;

    protected:
        /**
         * Generate a VL spec.
         * @return Generated VL spec.
         */
        vl::Spec spec() const override;

    private:
        /**
         * A local copy of the Verlet list.
         */
        Pairs exclPairs;
    };
}
