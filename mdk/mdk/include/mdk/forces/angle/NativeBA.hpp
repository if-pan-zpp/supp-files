#pragma once
#include "../../data/Primitives.hpp"
#include "../../simul/SimulVar.hpp"
#include "../../utils/Units.hpp"

namespace mdk {
    /**
     * "Native" part of the bond angle potential, applied when a triple has
     * an associated native bond angle.
     */
    class NativeBA: public SimulVar {
    private:
        friend class BondAngles;
        /**
         * isNative[i] is 1 if (i-1, i, i+1) has an associated native bond angle,
         * 0 otherwise.
         */
        Bytes isNative;

        /**
         * theta0[i] contains the native bond angle whenever such is defined.
         * For triples (i-1, i, i+1) which have no such associated angle, the
         * value is undefined.
         */
        Scalars theta0;

        /// Value of $k_\theta$.
        double CBA = 30.0 * eps/pow(rad, 2);

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
            auto diff = theta - theta0[i];
            V += CBA * diff * diff;
            dV_dth += 2.0 * CBA * diff;
        }
    };
}