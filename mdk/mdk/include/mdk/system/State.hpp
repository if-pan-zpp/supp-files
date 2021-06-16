#pragma once
#include "../data/Primitives.hpp"
#include "../utils/Topology.hpp"
#include "../model/Model.hpp"
#include "../simul/SimulVar.hpp"

namespace mdk {
    /**
     * An object containing the dynamical state of the simulation,
     * i.e. forces and the potential energy.
     */
    struct Dynamics {
        double V = 0.0;
        Vectors F;

        void zero(int n) {
            V = 0.0;
            F = Vectors(n, Vector::Zero());
        }
    };

    /**
     * An object containing the physical state of the simulation, i.e.
     * current positions, velocities of the residues, time and the box shape.
     */
    class State: public SimulVar {
    public:
        int n;
        Vectors r, v;
        double t;
        Topology top;

        Dynamics dyn;

        void prepareDyn();
        void updateWithDyn(Dynamics const& othDyn);
        void bind(Simulation& simul) override;

        /**
         * Set the positions and velocities of the residues in a \p Model
         * object to the current values.
         * @param model
         */
        void exportTo(Model& model) const;
    };
}
