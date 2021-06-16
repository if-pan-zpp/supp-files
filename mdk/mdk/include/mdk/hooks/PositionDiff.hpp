#pragma once
#include "Hook.hpp"
#include "../simul/Simulation.hpp"
#include "../system/State.hpp"
#include "../model/Model.hpp"
#include <fstream>
#include <string>
#include <map>

namespace mdk {
    /**
     * A hook which saves position diffs. I _assume_ this takes the path to the
     * positions as outputted by the Fortran implementation and compares the
     * positions in the consecutive steps with the current positions, whereafter
     * it outputs the diffs to the \p outputPath.
     */
    class PositionDiff : public Hook, SimulVar {
    public:
        /**
         * Create a \p PositionDiff object with given input and output paths.
         * @param inputPath Path to the file with Fortran positions.
         * @param outputPath Path to the diff file to create.
         */
        PositionDiff(std::string inputPath, std::string const& outputPath):
            inputPath(inputPath), output(outputPath, std::ofstream::out) {};

        /**
         * Bind the hook to the simulation; also fetch the state and load the
         * reference positions.
         * @param simulation
         */
        void bind(Simulation& simulation) override;

        /**
         * Execute the hook, i.e. output the diffs for the current step to the
         * diff file.
         * @param step_nr Number of the step of the simulation at which the hook
         * was invoked.
         */
        void execute(int step_nr) override;
    private:
        /**
         * A path to the Fortran reference positions file.
         */
        std::string inputPath;

        /**
         * Position diffs output stream.
         */
        std::ofstream output;

        /**
         * Reference positions loaded from the file at \p inputPath, in the
         * form of a map from the step number to a list of positions.
         */
        std::map<int, Vectors> refPositions;

        /**
         * Saved const pointer to the state of the simulation.
         */
        State const* state = nullptr;
    };
}
