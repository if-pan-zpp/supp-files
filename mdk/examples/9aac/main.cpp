#include <mdk/files/seq/LegacyParser.hpp>
#include <mdk/files/pdb/Parser.hpp>
#include <fstream>
using namespace mdk;
using namespace std;

int main() {
    auto model = seq::LegacyParser().read("data/9AAC.txt").asModel();

    auto rand = Random(0);
    model.morphIntoSAW(rand, false);

    ofstream coarseFile("data/9AAC.pdb");
    pdb::Parser().write(coarseFile, pdb::Data(pdb::Model(model)));
    return 0;
}
