#pragma once
#include "Force.hpp"
#include "../verlet/List.hpp"
#include "../verlet/Spec.hpp"

namespace mdk {
    /**
     * A nonlocal force interface. It differs from the standard force
     * interface in that it must provide some specs for the VL list, and
     * needs to implement a hook to run when a VL list is updated.
     */
    class NonlocalForce: public Force {
    protected:
        /**
         * A copy of the generated \p spec(), for use in the computation of
         * forces (for example it holds the cutoff distance).
         */
        mutable vl::Spec savedSpec;

        /**
         * Compute the spec to register to the VL list.
         * @return Generated VL spec.
         */
        virtual vl::Spec spec() const = 0;

        /**
         * A pointer to the VL list. It's not const chiefly because the update
         * hooks may modify the list, in particular create exclusions, like
         * for example the native contacts' forcefield removing the native
         * contacts
         */
        vl::List *vl = nullptr;

        /**
         * Registers the spec into the Verlet list.
         * Note: it must be executed only one all the data required for its
         * generation has been fetched, i.e. usually at the end of \p bind.
         */
        void installIntoVL();

    public:
        /**
         * Bind the nonlocal force to the simulation object. This class
         * shouldn't be bound to the simulation object, rather invoked by
         * derived classes. In particular, it fetches Verlet list into \p vl.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * An action to be performed when Verlet list is updated; usually one
         * copies the relevant pairs into a local list, perhaps with some
         * modifications. The actions are invoked in the order of being
         * added to the Verlet list.
         */
        virtual void vlUpdateHook() = 0;
    };
}