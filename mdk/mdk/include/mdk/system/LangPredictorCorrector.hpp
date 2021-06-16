#pragma once
#include "Integrator.hpp"
#include "../utils/Units.hpp"
#include "../data/Masses.hpp"

namespace mdk {
    /**
     * Langevin predictor-corrector fifth-order integrator, combined with a
     * Langevin noise. Apparently these two are combined because in the
     * Fortran code, the velocities are _set_ in \p lang or \p lang_mass. Now,
     * this looks like a bug, but we replicate it nonetheless.
     */
    class LangPredictorCorrector: public Integrator {
    public:
        explicit LangPredictorCorrector(double dt):
            dt(dt) {}

        void bind(Simulation& simulation) override;
        void init() override;
        void integrate() override;

        /**
         * Value of gamma for the Langevin noise.
         */
        double gamma = 2.0 * f77mass / tau;

    private:
        double dt;
        Masses m;
        // yi is 1/i! d^i r/dt^i from what I recall
        Vectors y0, y1, y2, y3, y4, y5;

        Random *random = nullptr;
        std::vector<Random> rngs;
        Vectors gaussianNoise;
        void generateNoise();

        bool initialized = false;

        double temperature = 0.35 * eps_kB;
    };
}
