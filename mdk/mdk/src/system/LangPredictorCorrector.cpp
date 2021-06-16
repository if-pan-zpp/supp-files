#include "system/LangPredictorCorrector.hpp"
#include "system/State.hpp"
#include "simul/Simulation.hpp"
using namespace mdk;

void LangPredictorCorrector::init() {
    for (int i = 0; i < state->n; ++i) {
        y2[i] = state->dyn.F[i]/m[i] * (dt*dt/2.0);
    }
    initialized = true;
}

void LangPredictorCorrector::generateNoise() {
    if (initialized) {
        #ifdef LEGACY_MODE
            #pragma omp task
            for (int dim = 0; dim < 3; ++dim) {
                for (int i = 0; i < state->n; ++i) {
                    gaussianNoise[i](dim) = random -> normal();
                }
            }
        #else
            for (int dim = 0; dim < 3; ++dim) {
                #pragma omp task
                {
                    for (int i = 0; i + 1 < state->n; i += 2) {
                        std::pair<double, double> normals = rngs[dim].two_normals();
                        gaussianNoise[i](dim) = normals.first;
                        gaussianNoise[i + 1](dim) = normals.first;
                    }
                    if (state->n % 2) {
                        gaussianNoise[state->n-1](dim) = rngs[dim].normal();
                    }
                }
            }
        #endif
    }
}

void LangPredictorCorrector::integrate() {
    double noiseVariance = sqrt(2.0*temperature *gamma*dt) * dt;
    double gamma_dt = gamma / dt;

    #pragma omp parallel for
    for (int i = 0; i < state->n; ++i) {
        // Damping and white noise
        y1[i] += gaussianNoise[i] * noiseVariance / m[i];
        state->dyn.F[i] -= gamma_dt * y1[i];

        // Correct
        Vector err = y2[i] - state->dyn.F[i]/m[i] * (dt*dt/2.0);
        y0[i] -= 3.0/16.0 * err;
        y1[i] -= 251.0/360.0 * err;
        y2[i] -= 1.0 * err;
        y3[i] -= 11.0/18.0 * err;
        y4[i] -= 1.0/6.0 * err;
        y5[i] -= 1.0/60.0 * err;


        // Predict
        y0[i] += y1[i] + y2[i] + y3[i] + y4[i] + y5[i];
        y1[i] += 2.0*y2[i] + 3.0*y3[i] + 4.0*y4[i] + 5.0*y5[i];
        y2[i] += 3.0*y3[i] + 6.0*y4[i] + 10.0*y5[i];
        y3[i] += 4.0*y4[i] + 10.0*y5[i];
        y4[i] += 5.0*y5[i];


        state->r[i] = y0[i];
        state->v[i] = y1[i]/dt;
    }

    state->t += dt;
}

void LangPredictorCorrector::bind(Simulation &simulation) {
    Integrator::bind(simulation);

    m = simulation.data<Masses>();
    random = &simulation.var<Random>();
    auto model = simulation.data<Model>();

    y0 = y1 = y2 = y3 = y4 = y5 = Vectors(model.n, Vector::Zero());
    for (int i = 0; i < model.n; ++i) {
        auto& res = model.residues[i];
        y0[i] = res.r;
        y1[i] = res.v * dt;
    }

    gaussianNoise = Vectors(model.n);
    simulation.addAsyncTask([this]() { this->generateNoise(); });

    #ifndef LEGACY_MODE
    Random r(*random);
    for (int i = 0; i < 3; ++i) {
        r = r.getNewRandom();
        rngs.push_back(r);
    }
    #endif
}
