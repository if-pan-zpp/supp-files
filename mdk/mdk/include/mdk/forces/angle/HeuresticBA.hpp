#pragma once
#include "../../data/Primitives.hpp"
#include "../../simul/SimulVar.hpp"
#include "../../utils/PairType.hpp"

namespace mdk {
    /**
     * "Heurestic" part of the bond angle potential, applied when a triple has
     * no associated native bond angle.
     */
    class HeuresticBA: public SimulVar {
    private:
        friend class BondAngles;
        /// Degree of the polynomial.
        static constexpr const int D = 6;

        /// Polynomial coefficients.
        double coeff[numOfPTs][D+1];

        /**
         * List of "types of triples", or rather pair types for (i, i+1) for
         * a triple (i-1, i, i+1).
         */
        Bytes angleTypes;

    public:
        /**
         * Bind the part to the simulation - in particular, find (and possibly
         * create) \p BondAngles force field and add itself to it.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Term of the formula for bond angle potential.
         * Note: it's inline in order for the compiler to inline it in
         * \p BondAngles.cpp file.
         * @param i The index i in the triple (i-1, i, i+1).
         * @param theta Value of the bond angle between i-1, i and i+1.
         * @param V Potential energy reference to add to.
         * @param dV_dth Derivative of potential energy wrt the angle theta
         * to add to.
         */
        void term(int i, double theta, double& V, double& dV_dth) const {
            double const* coeffs = coeff[angleTypes[i]];
            double V_loc = 0.0;
            double dV_dth_loc = 0.0;
            // Here we compute the polynomial with Horner scheme.
            for (int d = D; d >= 0; --d) {
                if (d > 0) dV_dth_loc = d * coeffs[d] + theta * dV_dth_loc;
                V_loc = coeffs[d] + theta * V_loc;
            }

            V += V_loc;
            dV_dth += dV_dth_loc;
        }
    };
}
