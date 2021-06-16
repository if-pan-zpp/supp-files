#include "hooks/ExportPDB.hpp"
#include "simul/Simulation.hpp"
#include "files/pdb/Model.hpp"
#include "files/pdb/Parser.hpp"
#include "files/pdb/Record.hpp"
#include <fstream>
#include <sstream>

using namespace mdk;
using namespace mdk::pdb;

void mdk::ExportPDB::bind(Simulation &simulation) {
    base = simulation.data<Model>();
    state = &simulation.var<State>();
}

void ExportPDB::execute(int step_nr) {
    if (state->t - tprev >= period) {
        auto remark = records::Remark();
        std::stringstream ss;
        ss << "TIME = " << state->t << " V = " << state->dyn.V << '\n';
        remark.number = 6;
        remark.text = ss.str();
        data << remark;

        auto modelRec = records::Model();
        modelRec.serialNum = modelIdx++;
        data << modelRec;

        state->exportTo(base);
        auto atoms = pdb::Data::onlyAtoms(pdb::Model(base));
        data << atoms;

        data << records::Endmdl();

        tprev = state->t;
    }
}

ExportPDB::~ExportPDB() {
    data << records::End();
    auto modelFile = std::ofstream(modelPath);
    pdb::Parser().write(modelFile, data);
}
