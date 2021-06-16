#include "system/Integrator.hpp"
using namespace mdk;

void Integrator::bind(Simulation &simulation) {
    state = &simulation.var<State>();
}
