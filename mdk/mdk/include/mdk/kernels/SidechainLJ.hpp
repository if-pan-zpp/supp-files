#pragma once
#include "../data/Primitives.hpp"
#include "LennardJones.hpp"

namespace mdk {
    /**
     * Sidechain version of the L-J potential (for a single pair of
     * \p depth and \p sink_max). In the CPC14.pdf it's described as an L-J
     * potential with a well bounded from both ends; in practice this is
     * implemented as a one-sided well, i.e. for r < sink_max, V = - depth and
     * the other side of the well is handled by the Pauli exclusion.
     */
    class SidechainLJ {
    public:
        double depth = 1.0 * eps;
        double sink_max = 10.0 * angstrom;

        SidechainLJ() = default;
        SidechainLJ(double depth, double sink_max):
            depth(depth), sink_max(sink_max) {};

        inline double cutoff() const {
            return LennardJones(sink_max, depth).cutoff();
        }

        /**
         * Compute the potential energy of the force field.
         * @param norm Distance between the residues.
         * @param V Variable to add the potential to.
         * @param dV_dn Variable to add the derivative to.
         */
        inline void computeV(double norm, double& V, double& dV_dn) const {
            if (norm <= sink_max) V -= depth;
            else LennardJones(sink_max, depth).computeV(norm, V, dV_dn);
        }

        /**
         * Compute and add the L-J force between two residues. The
         * templates are here in order for us to be able to pass Eigen
         * expressions to it.
         * @tparam T1 Type of an lvalue to add the force on the first residue
         * to.
         * @tparam T2 Type of an lvalue to add the force on the second residue
         * to.
         * @param unit Normalized vector between the residues.
         * @param norm Distance between the residues.
         * @param V Variable to add the potential to.
         * @param F1 Lvalue to add the force on the first residue to.
         * @param F2 Lvalue to add the force on the second residue to.
         */
        template<typename T1, typename T2>
        inline void computeF(VRef unit, double norm, double& V,
            T1 F1, T2 F2) const {

            double dV_dn = 0.0;
            computeV(norm, V, dV_dn);
            F1 -= dV_dn * unit;
            F2 += dV_dn * unit;
        }
    };
}
