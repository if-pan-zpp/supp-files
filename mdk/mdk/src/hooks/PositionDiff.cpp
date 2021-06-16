#include "hooks/PositionDiff.hpp"
#include "mdk/utils/Units.hpp"
#include <iostream>
#include <iomanip>
using namespace mdk;
using namespace std;

void PositionDiff::bind(Simulation& simulation) {
    state = &simulation.var<State>();

    auto& model = simulation.data<Model>();
    string token;
    std::ifstream input(inputPath, std::ifstream::in);

    while (input >> token) {
        if (token != "STEP") continue;
        int step_nr;
        input >> step_nr;

        Vectors r(model.n);

        for (int i = 0; i < model.n; ++i) {
            input >> r[i][0] >> r[i][1] >> r[i][2];
        }

        refPositions[step_nr] = r;
    }
}

void PositionDiff::execute(int step_nr) {
    auto it = refPositions.find(step_nr);
    if (it != refPositions.end()) {
        VectorBase diffs = (it -> second) - (state -> r / angstrom);
        Scalars diff_lengths = diffs.colwise().norm();
        double max_diff_len = diff_lengths.maxCoeff();
        double mean_diff_len = diff_lengths.mean();

        output << "Step #" << setw(8) << step_nr
               << ": max " << setw(12) << max_diff_len
               << ": mean " << setw(12) << mean_diff_len << endl;
    }
    else {
        output << "No information about positions in step " << step_nr << endl;
    }
}
