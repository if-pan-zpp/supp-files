#include "forces/dihedral/ComplexNativeDihedral.hpp"
#include "forces/dihedral/DihedralAngles.hpp"
using namespace mdk;

void ComplexNativeDihedral::bind(Simulation &simulation) {
    NativeDihedralBase::bind(simulation);
    auto& unifiedDih = simulation.var<DihedralAngles>();
    unifiedDih.natDih = this;
}