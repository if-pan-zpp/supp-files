#pragma once

namespace mdk {
    class Simulation;

    /**
     * A "simulation variable", i.e. an object which can be bound to a
     * particular simulation (usually to finish instantiation).
     */
    class SimulVar {
    public:
        /**
         * Bind an object to a simulation.
         * @param simulation Simulation to bind the object to.
         */
        virtual void bind(Simulation& simulation) = 0;
    };
}