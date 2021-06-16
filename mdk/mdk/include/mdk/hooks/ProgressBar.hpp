#pragma once
#include "../system/State.hpp"
#include "../simul/SimulVar.hpp"
#include "Hook.hpp"
#include <chrono>

namespace mdk {
    /**
     * A progress bar. Chiefly for convenience. Prints a status bar, current
     * (internal and clock) time and total simulation span, and the current
     * potential energy.
     * Note: Because it flushes the output, care must be taken not to invoke
     * it too often.
     */
    class ProgressBar: public Hook, SimulVar {
    public:
        /**
         * Construct a \p ProgressBar object.
         * @param totalTime Total span of the simulation; used for normalizing
         * the current time to display it as a percentage.
         * @param updatePeriod How often to update the progrss bar.
         * @param width Width of the progress bar in characters.
         */
        explicit ProgressBar(double totalTime,
            double updatePeriod, int width = 70);

        void bind(Simulation& simulation) override;

        void execute(int step_nr) override;

    private:
        using time_point = std::chrono::high_resolution_clock::time_point;

        Simulation *simul = nullptr;
        State *state = nullptr;

        /**
         * Real (clock) time when the simulation was started.
         */
        time_point realTime0;

        double totalTime, updatePeriod, prevTime;
        int width;
    };
}
