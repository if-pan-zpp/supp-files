#pragma once
#include "NativeDihedralBase.hpp"

namespace mdk {
    /**
     * The complex variant of the native part of the dihedral angle potential.
     */
    class ComplexNativeDihedral: public NativeDihedralBase {
    public:
        friend class DihedralAngles;
        double CDA = 0.66 * eps / pow(rad, 2.0);
        double CDB = 0.66 * eps / pow(rad, 2.0);

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
            double _phi0 = phi0[i];

            V += CDA * (1.0 - cos(phi - _phi0)) +
                 CDB * (1.0 - cos(3.0 * (phi - _phi0)));

            dV_dphi += CDA * sin(phi - _phi0) +
                       3.0 * CDB * sin(3.0 * (phi - _phi0));
        }
    };
}