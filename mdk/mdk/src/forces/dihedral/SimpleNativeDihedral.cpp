#include "forces/dihedral/SimpleNativeDihedral.hpp"
#include "forces/dihedral/DihedralAngles.hpp"
using namespace mdk;

void SimpleNativeDihedral::bind(Simulation &simulation) {
    NativeDihedralBase::bind(simulation);
    auto& unifiedDih = simulation.var<DihedralAngles>();
    unifiedDih.natDih = this;
}