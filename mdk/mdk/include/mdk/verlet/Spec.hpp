#pragma once

namespace mdk::vl {
    /**
     * A specification (spec for short) for the non-local force. Technically
     * one could generate a Verlet list for every non-local force, however this
     * would be wasteful - thus we generate a most restrictive possible Verlet
     * list which contains the "sub-Verlet lists" of the particular non-local
     * forces and let them filter the general list on their own accord into a
     * local list.
     */
    struct Spec {
        /**
         * A square of the cutoff distance.
         */
        double cutoffSq = 0.0;

        /**
         * Minimum bond separation between the residues of one pair in the
         * Verlet list. Technically this is just an optimization.
         */
        int minBondSep = 3;
    };
}
