#pragma once
#include "../../simul/SimulVar.hpp"
#include "../../simul/Simulation.hpp"

namespace mdk {
    /**
     * A common part of simple and complex native dihedral variants. In
     * particular, it stores whether a quadruple (i-2, i-1, i, i+1) has an
     * associated native dihedral angle, and its value if it is the case.
     */
    class NativeDihedralBase: public SimulVar {
    protected:
        /**
         * isNative[i] = 1 if a quadruple (i-2, i-1, i, i+1) has an associated
         * native dihedral angle, 0 otherwise.
         */
        Bytes isNative;

        /**
         * phi0[i] has the value of an associated native dihedral angle, if such
         * exists for the quadruple (i-2, i-1, i, i+1); otherwise, the value is
         * undefined.
         */
        Scalars phi0;

    public:
        /**
         * Bind the variant to the simulation. Technically one shouldn't
         * add this base to the simulation itself; it serves as a prototype for
         * \p bind in \p ComplexNativeDihedral and \p SimpleNativeDihedral --
         * in particular, it retrieves and computes \p isNative and \p phi0.
         * @param simulation
         */
        void bind(Simulation& simulation) override;
    };
}