#pragma once
#include "../utils/Units.hpp"
#include "../system/State.hpp"
#include "../model/Model.hpp"
#include "../simul/SimulVar.hpp"
#include "../files/pdb/Data.hpp"
#include "Hook.hpp"
#include <filesystem>
#include <limits>

namespace mdk {
    /**
     * The hook for exporting current positions of the residues to
     * a PDB file. The positions are stored as MODELs within the PDB file;
     * (at this moment) the file is output only after the resolution of the
     * simulation, i.e. in the destructor. Also, at the moment only the
     * positions (and not for example REMARKs with QA stats or contacts) are
     * saved.
     */
    class ExportPDB: public Hook, SimulVar {
    private:
        /**
         * Reference model to export the current positions to.
         */
        Model base;

        /**
         * A PDB data object containing the currently stored records (in
         * particular the positions).
         */
        pdb::Data data;

        /**
         * Current model serial number (1-indexed).
         */
        int modelIdx = 1;

        /**
         * Const pointer to the state.
         */
        State const* state = nullptr;

        /**
         * Path to the PDB file.
         */
        std::filesystem::path modelPath;

        /**
         * Span betwen consecutive exports (in internal time).
         */
        double period;

        /**
         * Time point when the record was last exported (in internal time).
         * This one is initialized to minimum value so that the first
         * execution always exports.
         */
        double tprev = std::numeric_limits<double>::lowest();

    public:
        /**
         * Construct a \p ExportPDB hook with a model path and the period.
         * @param modelPath Path to the file to which to output the positions.
         * @param period Span between consecutive exports.
         */
        explicit ExportPDB(std::filesystem::path modelPath,
            double period = 1000.0 * microsecond):
            modelPath(std::move(modelPath)),
            period(period) {};

        /**
         * This destructor saves the \p to the file.
         */
        ~ExportPDB();

        /**
         * Bind the hook to the simulation. Also, fetch the state.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Execute the hook, i.e. check if the span since the last time the
         * positions were exported is greater than \p period and if it is,
         * export the positions and set the last time to current time.
         * @param step_nr Number of the step of the simulation when the
         * hook is invoked.
         */
        void execute(int step_nr) override;
    };
}
