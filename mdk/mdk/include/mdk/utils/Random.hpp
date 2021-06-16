#pragma once
#include <vector>
#include <random>
#include <Eigen/Core>

namespace mdk {
    /**
     * A random number generator. We use two versions: a legacy version taken
     * from Fortran, and a modern version. Aside from sampling from [0, 1], we
     * add sampling from [a, b], N(0, 1), N(mu, sigma^2) or points on S^2.
     */
    class Random {
#ifdef LEGACY_MODE
    private:
        static constexpr int
            im1 = 2147483563,
            im2 = 2147483399,
            imm1 = im1-1,
            ia1 = 40014,
            ia2 = 40692,
            iq1 = 53668,
            iq2 = 52774,
            ir1 = 12211,
            ir2 = 3791,
            ntab = 32,
            ndiv = 1+imm1/ntab;

        static constexpr double
            eps = 1.2e-7,
            rnmx = 1.-eps,
            am = (float) 1.0/im1;

        int iy = 0, idum = -448, idum2 = 123456789;
        int iv[ntab];
    public:
        Random(int seed) {
            idum = -seed;
            for (auto& x: iv) x = 0;
        }
        inline double uniform() {
            int k, j;
            if (idum <= 0) {
                idum2 = idum = std::max(-idum, 1);
                for (j = ntab + 7; j >= 0; --j) {
                    k = idum / iq1;
                    idum = ia1 * (idum - k * iq1) - k * ir1;
                    if (idum < 0) idum += im1;
                    if (j < ntab) iv[j] = idum;
                }
                iy = iv[0];
            }

            k = idum / iq1;
            idum = ia1 * (idum - k * iq1) - k * ir1;
            if (idum < 0) idum += im1;

            k = idum2 / iq2;
            idum2 = ia2 * (idum2 - k * iq2) - k * ir2;
            if (idum2 < 0) idum2 += im2;

            j = iy / ndiv;
            iy = iv[j] - idum2;
            iv[j] = idum;
            if (iy < 1) iy += imm1;

            return std::min(am * iy, rnmx);
        }
#else
    private:
        uint64_t state;
        Random() = default;

    public:
        Random(int seed) {
            assert (seed > 0);
            state = seed;

            // Shuffle the state a bit
            for (int i = 0; i < 100; ++i) uniform();
        }

        inline double uniform() {
            static const double inv = 1.0 / (double) (1ull << 32);
            uint64_t result = state * 0xd989bcacc137dcd5ull;
            state ^= state >> 11;
            state ^= state << 31;
            state ^= state >> 18;
            double res = (uint32_t) (result >> 32ull);
            return res * inv;
        }

        /* This function is used to create a new random object
           that won't generate similar values.
           Uses splitmix64 to create new state */
        Random getNewRandom() {
            uint64_t result = state;
            result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
            result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
            result = result ^ (result >> 31);
            Random rng;
            rng.state = result;
            return rng;
        }

#endif
        Random(Random const& oth) = default;

        inline double uniform(double a, double b) {
            return a + (b - a) * uniform();
        }

        inline double normal() {
            double r1 = uniform();
            double r2 = uniform();
            return sqrt(-2.0  * log(r1)) * cos(2.0 * M_PI * r2);
        }

        inline std::pair<double, double> two_normals() {
            double r1 = uniform();
            double r2 = uniform();
            double r = sqrt(-2.0  * log(r1));
            double ang = 2.0 * M_PI * r2;
            return {r * cos(ang), r * sin(ang)};
        }

        inline double normal(double mu, double sigma) {
            return mu + sigma * normal();
        }

        inline Eigen::Vector3d sphere() {
            Eigen::Vector3d r { normal(), normal(), normal() };
            return r.normalized();
        }
    };
}
