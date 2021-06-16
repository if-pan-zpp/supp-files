#include <mdk/files/seq/LegacyParser.hpp>
#include <mdk/files/pdb/Parser.hpp>
#include <fstream>
#include <iostream>
using namespace mdk;
using namespace std;

int main() {
    auto model = seq::LegacyParser().read("data/glut.txt").asModel();

    auto rand = Random(0);
    model.morphIntoSAW(rand);
    cout << model.residues[1000].r << '\n';

    ofstream coarseFile("data/glut.pdb");
    auto data = pdb::Data::onlyAtoms(pdb::Model(model));
    pdb::Parser().write(coarseFile, data);
    return 0;
}
