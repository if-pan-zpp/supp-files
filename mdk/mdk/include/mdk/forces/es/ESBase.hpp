#pragma once
#include "../NonlocalForce.hpp"

namespace mdk {
    /**
     * A base class for the Debye-Hueckel screened electrostatic potentials.
     */
    class ESBase: public NonlocalForce {
    protected:
        /**
         * A struct detailing a contact between pairs. We store the indices
         * and the product of the charges.
         */
        struct Contact {
            /// Index of the first residue.
            int i1;
            /// Index of the second residue.
            int i2;
            /// Product of residues' charges.
            double q1_x_q2;
        };

        /**
         * A local, filtered Verlet list of pairs. Since most residues are not
         * charged, and pairs even moreso, it benefits us to keep a local list
         * and not check the pairs which we would know cannot interact with
         * each other.
         */
        std::vector<Contact> pairs;

        /**
         * A list of charges of the residues (in units of \p echarge).
         */
        Eigen::Matrix<int8_t, Eigen::Dynamic, 1> charge;

    public:
        /**
         * Bind the class to the simulation. It initializes \p charge and adds
         * the force to the VL list.
         * Note: this class (\p ESBase) shouldn't be added to the simulation
         * class, only the derived classes.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * An action to be performed when the Verlet list is updated. Here we
         * want to recreate \p pairs, i.e. filter the Verlet list so that only
         * charged pairs are present.
         */
        void vlUpdateHook() override;
    };
}
