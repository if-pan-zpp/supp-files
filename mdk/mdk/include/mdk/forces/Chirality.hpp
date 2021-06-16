#pragma once
#include "Force.hpp"

namespace mdk {
    /**
     * Chirality potential for the models with a native structure (in
     * particular the positions of the residues). It's supposed to be a more
     * replacement for bond and dihedral potentials, if the native structure
     * is provided.
     */
    class Chirality: public Force {
    private:
        /**
         * A list of cube inverses of $d_0$, where $d_0$ is the length of
         * $v_i = r_{i+1} - r_i$ in the native structure.
         */
        Scalars d0_cube_inv;

        /**
         * A list of the values of C_i in the native structure.
         */
        Scalars C_nat;

        /**
         * inRange[i] = 1 if a quadruple (i-2, i+1, i, i+1) is connected,
         * i.e. in a single chain.
         */
        Bytes inRange;

    public:
        /**
         * The amplitude of the potential.
         */
        double e_chi = 1.0 * eps;

        /**
         * Bind the force field to the simulation.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Asynchronous part of the computation.
         * @param dynamics Dynamics object to add potential energy and forces
         * to.
         */
        void asyncPart(Dynamics &dynamics) override;
    };
}
