#include <mdk/simul/Simulation.hpp>
#include <mdk/files/pdb/Parser.hpp>
#include <mdk/files/param/LegacyParser.hpp>
#include <mdk/system/LangPredictorCorrector.hpp>
#include <mdk/forces/All.hpp>
#include <mdk/hooks/ProgressBar.hpp>
#include <mdk/hooks/ExportPDB.hpp>
#include <fstream>
#include <omp.h>
using namespace mdk;
using namespace std;

int main() {
    omp_set_num_threads(1);

    ifstream pdb_file("data/1ubq.pdb");
    auto atomic = pdb::Parser().read(pdb_file).asModel();

    ifstream param_file("data/parametersMJ96.txt");
    auto params = param::LegacyParser().read(param_file);

    atomic.addContactsFromAtomOverlap();
    auto model = atomic.coarsen();

    auto rand = Random(448);
    rand.uniform();

    model.legacyMorphIntoSAW(rand, false, 0.0, 4.56*angstrom, true);
    model.initVelocity(rand, 0.35 * eps_kB, false);

    Simulation simul(model, params);

    simul.add<Random>(rand);
    simul.add<LangPredictorCorrector>(0.005 * tau);

    simul.add<Tether>(true);
    simul.add<NativeBA>();
    simul.add<ComplexNativeDihedral>();

    simul.add<NativeContacts>();
    simul.add<PauliExclusion>();

    auto total = 15000.0*tau;
    simul.add<ProgressBar>(total, 10.0*tau);
    simul.add<ExportPDB>("model.x.pdb", 100.0*tau);

    simul.step(total);

    return 0;
}
