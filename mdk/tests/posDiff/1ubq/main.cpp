#include <mdk/simul/Simulation.hpp>
#include <mdk/files/pdb/Parser.hpp>
#include <mdk/files/param/LegacyParser.hpp>
#include <mdk/hooks/PositionDiff.hpp>
#include <mdk/system/Leapfrog.hpp>
#include <mdk/system/LangPredictorCorrector.hpp>
#include <mdk/forces/Tether.hpp>
#include <mdk/forces/angle/NativeBA.hpp>
#include <mdk/forces/dihedral/ComplexNativeDihedral.hpp>
#include <mdk/forces/go/NativeContacts.hpp>
#include <mdk/forces/PauliExclusion.hpp>
#include <fstream>
using namespace mdk;
using namespace std;

int main() {
    ifstream pdb_file("data/1ubq.pdb");
    auto atomic = pdb::Parser().read(pdb_file).asModel();

    ifstream param_file("data/parametersMJ96.txt");
    auto params = param::LegacyParser().read(param_file);

    atomic.addContactsFromAtomOverlap();
    auto model = atomic.coarsen();

    auto rand = Random(448);
    rand.uniform();

    model.legacyMorphIntoSAW(rand, false, 0, 4.56*angstrom, true);
    model.initVelocity(rand, 0.35 * eps_kB, false);

    Simulation simul(model, params);

    simul.add<Random>(rand);
    simul.add<LangPredictorCorrector>(0.005 * tau);

    simul.add<PositionDiff>("data/positions.txt", "posDiffs.out");

    simul.add<Tether>(true);
    simul.add<NativeBA>();
    simul.add<ComplexNativeDihedral>();

    simul.add<NativeContacts>();
    simul.add<PauliExclusion>();

    for (int i = 0; i < 100; ++i) {
    	simul.step();
    }
    return 0;
}
