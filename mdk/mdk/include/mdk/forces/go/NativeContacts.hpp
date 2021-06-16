#pragma once
#include "../NonlocalForce.hpp"

namespace mdk {
    /**
     * Go model potential. It is somewhat special as it interferes with
     * the Verlet list and with the quasi-adiabatic potential. It should go
     * first in the list of nonlocal forces added to the simulation object.
     */
    class NativeContacts: public NonlocalForce {
    public:
        /**
         * A struct with contact data.
         */
        struct Contact {
            /// First residue.
            int i1;

            /// Second residue.
            int i2;

            /// Lennard-Jones potential minimal distance.
            double r_min;
        };

        /**
         * Bind the force to the simulation, and also to the Verlet list.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation &simulation) override;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        void asyncPart(Dynamics &dynamics) override;

        /**
         * Action to perform when the Verlet list is updated. We want to (a)
         * retrieve the pairs that are in native contact into a list (\p
         * curPairs), (b) remove those pairs from the global Verlet list.
         */
        void vlUpdateHook() override;

    private:
        /**
         * Generate a spec for the Verlet list.
         * @return Generated spec.
         */
        vl::Spec spec() const override;

        /**
         * A "copy" of the global Verlet list to be swapped with it when we
         * update it; we want to do the swap so as to avoid excessive allocation
         * of memory when the program runs.
         */
        Pairs newVL;

        /// List of all native contacts.
        std::vector<Contact> allContacts;

        /// List of native contacts that are within the cutoff distance.
        std::vector<Contact> curPairs;

        /**
         * Depth of the Lennard-Jones potential, with which the natively
         * connected residues interact.
         */
        double depth = 1.0 * eps;
    };
}
