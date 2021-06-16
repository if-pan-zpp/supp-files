#pragma once
#include "NativeDihedralBase.hpp"

namespace mdk {
    /**
     * The simple variant of the native part of the dihedral angle potential.
     */
    class SimpleNativeDihedral: public NativeDihedralBase {
    public:
        friend class DihedralAngles;
        double CDH = 3.33 * eps/pow(rad, 2);

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
            auto diff = phi - phi0[i];
            V += 0.5 * CDH * diff * diff;
            dV_dphi += CDH * diff;
        }
    };
}