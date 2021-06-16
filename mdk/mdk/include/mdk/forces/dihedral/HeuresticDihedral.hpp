#pragma once
#include "../../files/param/Parameters.hpp"
#include "../../simul/SimulVar.hpp"
#include "../../simul/Simulation.hpp"

namespace mdk {
    /**
     * Heurestic part of the dihedral potential, applied when with the
     * quadruple (i-2, i-1, i, i+1) is not associated a native dihedral angle.
     */
    class HeuresticDihedral: public SimulVar {
    private:
        friend class DihedralAngles;
        /**
         * Coefficients per pair type (where pair is (i-1, i) from the
         * quadruple.
         */
        double coeff[numOfPTs][6];

        /**
         * Types of angles (i.e. pair types for each pair (i-1, i)).
         */
        Eigen::Matrix<int8_t, Eigen::Dynamic, 1> angleTypes;

    public:
        /**
         * Bind the part to the simulation - in particular, find (and possibly
         * create) \p DihedralAngles force field and add itself to it.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Term of the formula for the dihedral angle potential of a single
         * quadruple.
         * Note: it's inline in order for the compiler to inline it in
         * \p DihedralAngles.cpp file.
         * @param i The index i in the quadruple (i-2, i-1, i, i+1).
         * @param phi Dihedral angle between planes defined by sequences
         * i-2, i-1, i and i-1, i, i+1.
         * @param V Potential energy reference to add to.
         * @param dV_dphi Derivative of potential energy wrt the angle phi
         * to add to.
         */
        void term(int i, double phi, double& V, double& dV_dphi) const {
            double sin_phi = sin(phi), cos_phi = cos(phi);
            double sin_2_phi = sin_phi * sin_phi;
            double cos_2_phi = cos_phi * cos_phi;

            double const* C = coeff[angleTypes[i]];

            V += C[0]
               + C[1] * sin_phi
               + C[2] * cos_phi
               + C[3] * sin_2_phi
               + C[4] * cos_2_phi
               + C[5] * sin_phi * cos_phi;

            dV_dphi += C[1] * cos_phi
                     - C[2] * sin_phi
                     + 2.0 * (C[3] - C[4]) * sin_phi * cos_phi
                     + C[5] * (cos_2_phi - sin_2_phi);
        }
    };
}
