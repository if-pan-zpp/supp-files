#include "forces/NonlocalForce.hpp"
#include "simul/Simulation.hpp"
using namespace mdk;

void NonlocalForce::bind(Simulation &simulation) {
    Force::bind(simulation);
    vl = &simulation.var<vl::List>();
}

void NonlocalForce::installIntoVL() {
    savedSpec = spec();
    vl->registerNF(*this, savedSpec);
}