#include "simul/Simulation.hpp"
#include "forces/Force.hpp"
#include "forces/NonlocalForce.hpp"
#include "hooks/Hook.hpp"
#include "system/Integrator.hpp"
using namespace mdk;

extern Dynamics thread_dyn;
#pragma omp threadprivate(thread_dyn)
Dynamics thread_dyn;

void Simulation::calcForces() {
    state -> prepareDyn();
    verlet_list -> check();

    #pragma omp parallel
    {
        thread_dyn.zero(state -> n);

        #pragma omp master
        for (auto const& task : asyncTasks) {
            task();
        }
            
        for (auto* force: forces) {
            force->asyncPart(thread_dyn);
        }

        #pragma omp critical
        {
            state -> updateWithDyn(thread_dyn);
        }
    }

    for (auto* force: forces) {
        force->syncPart(state -> dyn);
    }
}

void Simulation::init() {
    state = &var<State>();
    verlet_list = &var<vl::List>();
    
    step_nr = 0;

    calcForces();
    integrator->init();

    for (auto* hook: hooks) {
        hook->execute(0);
    }
    
    initialized = true;
}

void Simulation::step() {
    if (not initialized) {
        init();
    }
    
    step_nr++;

    calcForces();
    integrator->integrate();

    for (auto* hook: hooks) {
        hook->execute(step_nr);
    }
}

void Simulation::step(double t) {
    auto& state = var<State>();
    auto t0 = state.t;
    while (state.t - t0 < t) step();
}
