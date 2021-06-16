#pragma once
#include "Stat.hpp"
#include "../simul/SimulVar.hpp"
#include <vector>
#include "../data/Types.hpp"
#include "../files/param/Parameters.hpp"

namespace mdk {
    /**
     * An object containing stats for the residues. Also initializes
     * these numbers in one place.
     */
    class Stats: public SimulVar {
    private:
        /**
         * Types of amino acids.
         */
        Types const *types = nullptr;

        /**
         * List of polarization values taken from the \p param::Parameters
         * object.
         */
        std::vector<param::Polarization> polarization;

    public:
        /**
         * List of stats for each residue.
         */
        std::vector<Stat> stats;

        void bind(Simulation& simulation) override;

        /**
         * Type of a contact (backbone-backbone etc.).
         */
        enum class Type: int8_t {
            BB, BS, SB, SS
        };

        /**
         * Generate diffs for if a contact between residues \p i1 and \p i2
         * of type \p type were created. It's inline so as for the compiler to
         * insert it and not have to call the function.
         * @param i1 First residue.
         * @param i2 Second residue.
         * @param type Type of a contact.
         * @param diffs An array of diffs; diffs[0] is the diff for the first
         * residue, diffs[1] is the diff for the second residue.
         */
        inline void creationDiffs(int i1, int i2, Type type, Stat diffs[2]) const {
            diffs[0] = diffs[1] = Stat();

            switch (type) {
            case Type::BB:
                --diffs[0].backbone;
                --diffs[1].backbone;
                break;
            case Type::BS:
                --diffs[0].backbone;
                --diffs[1].sidechain;
                break;
            case Type::SB:
                --diffs[0].sidechain;
                --diffs[1].backbone;
                break;
            case Type::SS:
                --diffs[0].sidechain;
                --diffs[1].sidechain;
            }

            if (type == Type::SS) {
                int idx[2] = {i1, i2};
                for (int k = 0; k < 2; ++k) {
                    switch (polarization[(int8_t)(*types)[idx[1-k]]]) {
                    case param::Polarization::POLAR:
                    case param::Polarization::POLAR_NEG:
                    case param::Polarization::POLAR_POS:
                        --diffs[k].polarSS;
                        break;

                    case param::Polarization::HYDROPHOBIC:
                        --diffs[k].hydrophobicSS;
                        break;

                    default:
                        break;
                    }
                }
            }
        }
    };
}