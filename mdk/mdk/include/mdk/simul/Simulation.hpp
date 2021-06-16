#pragma once
#include "../model/Model.hpp"
#include "../files/param/Parameters.hpp"
#include "../data/DataFactory.hpp"
#include "../verlet/List.hpp"
#include "SimulVar.hpp"
#include <typeindex>
#include <type_traits>
#include <stdexcept>
#include <functional>

namespace mdk {
    // These forward decls are necessary to avoid include loops, as all of
    // them include this file.
    class Force;
    class NonlocalForce;
    class Hook;
    class Integrator;

    /**
     * The main class of the library, responsible for:
     * - storing the simulation state;
     * - storing the forces and other associated objects, and (optionally)
     *   binding them to the simulation;
     * - running, or rather stepping through, the simulation.
     */
    class Simulation {
    public:
        /**
         * Initialize simulation from a \p Model object and the parameters.
         * @param model Model of the simulation.
         * @param params Parameters of the simulation.
         */
        Simulation(Model model, param::Parameters params):
            model(std::move(model)), params(std::move(params)),
            df(this->model, this->params) {};

        /**
         * Access the "static" data as created by the \p DataFactory. We access
         * the data by the type, but insofar as these objects represent functions
         * of the state and the simulation, this is (in my estimation)
         * appropriate.
         * Providing an invalid type will lead to linker error.
         * @tparam Data Type of static data to access.
         * @return Const reference to the data object.
         */
        template<typename Data>
        Data const& data() {
            return df.data<Data>();
        }

        /**
         * Accesses the variables associated with the simulation. The access is
         * by type, thus it's most suited to variables that are effectively
         * functions of the simulation (for example the physical state or the
         * Verlet list). If the referenced variable has not been yet added,
         * the function tries to emplace it (if it's default-constructible) --
         * such is, for example, the case with \p State or \p vl::List. Cannot be
         * run after simulation initialization. The reference is never invalidated.
         * @tparam Var Type of the variable to access.
         * @return Accessed variable. The function returns a non-const reference,
         * thus care must be taken as to modify the variables only on agreed-upon
         * points of the simulation (say, integration of forces).
         */
        template<typename Var>
        Var& var() {
            if (initialized) {
                throw std::runtime_error("Cannot access state vars after initialization");
            }

            auto idx = std::type_index(typeid(Var));
            if (vars.find(idx) == vars.end()) {
                // The variable has not been found, we try to construct it.
                if constexpr (std::is_default_constructible_v<Var>) {
                    return add<Var>();
                }
                else {
                    static_assert("Cannot default-construct Var");
                }
            }
            else {
                // The variable was added before, we simply access it.
                auto& _var = vars.at(idx);
                return *(Var*)_var.get();
            }
        }

        /**
         * Add a variable to the simulation. Depending on the type of the
         * variable, additional actions may be performed (such as binding to the
         * simulation or adding forces and/or nonlocal forces to appropriate lists).
         * Cannot be run after simulation initialization. The reference is
         * never invalidated.
         * @tparam T Type of the variable to add.
         * @tparam Args Types of arguments to use in construction.
         * @param args Values of arguments to use in construction.
         * @return Reference to a constructed and emplaced \p T(\p Args...).
         */
        template<typename T, typename... Args>
        T& add(Args&&... args) {
            if (initialized) {
                throw std::runtime_error("Cannot add state vars after initialization");
            }

            // Create the variables as \p std::shared_ptr<T>
            auto _var = std::make_shared<T>(std::forward<Args>(args)...);

            // Get a type index for T and add \p _var to vars; get an iterator
            // and a reference to T&
            auto idx = std::type_index(typeid(T));
            auto varIter = vars.emplace(idx, _var).first;
            auto& x = *(T*)varIter->second.get();

            if constexpr (std::is_base_of_v<SimulVar, T>) {
                // If a \p SimulVar, bind it to the simulation
                ((SimulVar&)x).bind(*this);
            }

            if constexpr (std::is_base_of_v<Force, T>) {
                // If a \p Force, add it to forces' list
                forces.push_back(&(Force&)x);
            }

            if constexpr (std::is_base_of_v<NonlocalForce, T>) {
                // If a \p NonlocalForce, add it to nonlocal forces' list
                nonlocalForces.push_back(&(NonlocalForce&)x);
            }

            if constexpr (std::is_base_of_v<Hook, T>) {
                // If a \p Hook, add it to hooks' list
                hooks.push_back(&(Hook&)x);
            }

            if constexpr (std::is_base_of_v<Integrator, T>) {
                // If an \p Integrator, replace the current integrator by it
                integrator = &(Integrator&)x;
            }

            // Return the reference
            return x;
        }

        /**
         * Add an "async task".  Cannot be run after simulation initialization.
         * @param f Function/lambda to add. (Seemingly the lambda cannot
         * capture the environment, i.e. must be trivial)
         */
        inline void addAsyncTask(std::function<void()> const& f) {
            if (initialized) {
                throw std::runtime_error("Cannot add async tasks after initialization");
            }

            asyncTasks.push_back(f);
        }

        /**
         * Initialize the simulation.
         */
        void init();

        /**
         * Step once through the simulation; the time step is as determined by
         * the used integrator (usually 5 ps).
         */
        void step();

        /**
         * Step a total of time \p t through the simulation.
         * @param t Time to advance.
         */
        void step(double t);

    private:
        /// Model of the simulation.
        Model model;

        /// Parameters of the simulation.
        param::Parameters params;

        /// Data factory.
        DataFactory df;

        /**
         * A heterogeneous collection of objects (stored as shared pointers to
         * void) indexed by their types.
         */
        std::unordered_map<std::type_index, std::shared_ptr<void>> vars;

        /// Force fields of the simulation
        std::vector<Force*> forces;

        /**
         * List of nonlocal forces of the simulation - the distinction is used
         * when the Verlet list is updated, so as to invoke appropriate update
         * functions.
         */
        std::vector<NonlocalForce*> nonlocalForces;

        /// Hooks of the simulation
        std::vector<Hook*> hooks;

        /// Async tasks.
        std::vector<std::function<void()>> asyncTasks;

        /// Pointer to the integrator used in the simulation.
        Integrator *integrator = nullptr;

        State *state = nullptr;
        vl::List *verlet_list = nullptr;

        /**
         * Whether simulation has yet been initialized; used for running the
         * init function only once.
         */
        bool initialized = false;

        /// Number of steps performed.
        int step_nr = 0;

        /// Internal function for invoking force fields.
        void calcForces();
    };
}
