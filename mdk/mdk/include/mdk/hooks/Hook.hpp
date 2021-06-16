#pragma once

namespace mdk {
    /**
     * This interface defines actions to be performed after the force
     * evaluation, i.e. a hook.
     */
    class Hook {
    public:
        /**
         * Execute the hook.
         * @param step_nr Number of the step of the simulation when the
         * hook is activated.
         */
        virtual void execute(int step_nr) = 0;
    };
}
