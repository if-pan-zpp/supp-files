#pragma once
#include "../system/State.hpp"
#include "../simul/SimulVar.hpp"
#include "../simul/Simulation.hpp"

namespace mdk {
    /**
     * An abstract interface for the force fields. Certain forces (QA in
     * particular) have a multi-stage evaluation. Thus we separated the
     * computation to an "asynchronous" and "synchronous" parts, where the
     * latter get executed in order of addition to the \p Simulation object.
     */
    class Force: public SimulVar {
    protected:
        /**
         * A const pointer to the \p State of simulation. The const-ness
         * is a declaration that the \p State is only modified by the
         * integrator, and thus \p asyncPart can be safely executed.
         */
        State const* state = nullptr;

    public:
        /**
         * Bind the force to the simulation. This class shouldn't be actually
         * added to the simulation, but rather serves as a prototype for actual
         * forces; in particular it saved \p state from the \p Simulation object.
         * @param simulation
         */
        void bind(Simulation& simulation) override;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        virtual void asyncPart(Dynamics& dynamics);

        /**
         * Synchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and forces to.
         */
        virtual void syncPart(Dynamics& dynamics);
    };
}
