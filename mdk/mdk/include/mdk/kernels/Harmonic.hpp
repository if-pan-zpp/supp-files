#pragma once
#include "../utils/Units.hpp"
#include "../data/Primitives.hpp"

namespace mdk {
    /**
     * Harmonic force "kernel", i.e. a separated class responsible only for
     * the computation of the formula. The definitions are inlined in order
     * for the compiler to inline them.
     */
    class Harmonic {
    public:
        /**
         * "Harmonic" part of the harmonic potential.
         */
        double H1 = 50.0 * eps/pow(angstrom, 2.0);

        /**
         * "Anharmonic" part of the harmonic potential. From what I have seen
         * it's pretty much always zero, perhaps it's a legacy/deprecated
         * feature.
         */
        double H2 = 0.0;

        Harmonic() = default;
        Harmonic(double H1, double H2):
            H1(H1), H2(H2) {};

        /**
         * Compute the potential energy of the force field.
         * @param dx Displacement, i.e. the difference between current length
         * and the equilibrium length
         * @param V Variable to add the potential to.
         * @param dV_dx Variable to add the derivative to.
         */
        inline void computeV(double dx, double& V, double& dV_dx) const {
            auto dx2 = dx*dx;
            V += dx2 * (H1 + H2 * dx2);
            dV_dx += dx * (2.0 * H1 + 4.0 * H2 * dx2);
        }

        /**
         * Compute and add the harmonic force between two residues. The
         * templates are here in order for us to be able to pass Eigen
         * expressions to it.
         * @tparam T1 Type of an lvalue to add the force on the first residue
         * to.
         * @tparam T2 Type of an lvalue to add the force on the second residue
         * to.
         * @param unit Normalized vector between the residues.
         * @param dx Displacement, i.e. the difference between current length
         * and the equilibrium length;
         * @param V Variable to add the potential to.
         * @param F1 Lvalue to add the force on the first residue to.
         * @param F2 Lvalue to add the force on the second residue to.
         */
        template<typename T1, typename T2>
        inline void computeF(VRef unit, double dx, double& V,
            T1 F1, T2 F2) const {

            double dV_dn = 0.0;
            computeV(dx, V, dV_dn);
            F1 += dV_dn * unit;
            F2 -= dV_dn * unit;
        }
    };
}
