#pragma once
#include <vector>
#include <unordered_map>
#include <optional>
#include "../../utils/AminoAcid.hpp"
#include "../../utils/TupleHash.hpp"
#include "../../utils/PairType.hpp"

namespace mdk::param {
    /**
     * Polarization values, as derived from CPC14.pdf and from reading
     * cg.f.
     */
    enum class Polarization: int8_t {
        GLY_PRO = 0, HYDROPHOBIC = 1, MISSING = 2,
        POLAR = 3, POLAR_NEG = 4, POLAR_POS = 5
    };

    /** A POD representing the parameters from the parameter file. Arguably
     * the legacy file could be split so that everything wouldn't be in one
     * structure.
     */
    class Parameters {
    public:
        Parameters();

        using Coeffs = std::vector<double>;
        /**
         * "Default" angle parameters. These are I assume for the simplified
         * angle potential (option \p lsimpang).
         */
        Coeffs defAngleParams;

        template<typename T>
        using PerPTData = std::unordered_map<PairType, T>;

        /**
         * A set of bond angle parameters (for heurestic potential) for each
         * pair of types of residues i, i+1 in a triple (i-1, i, i+1).
         */
        PerPTData<Coeffs> angleParams;

        /**
         * A set of dihedral angle parameters (for heurestic potential) for each
         * pair of types of residues i-1, i in a quadruple (i-2, i-1, i, i+1).
         */
        PerPTData<Coeffs> dihedralParams;

        /// Specificity parameters POD.
        struct SpecificityParams {
            /// Polarization of an amino acid.
            Polarization polarization;

            /// Maximum number of sidechain contacts an amino acid can form.
            int maxSidechain;

            /**
             * Maximum number of sidechain contacts an amino acid can form with
             * hydrophobic amino acids.
             */
            int maxHydrophobicSS;

            /**
             * Maximum number of sidechain contacts an amino acid can form with
             * polar amino acids.
             */
            int maxPolarSS;
        };

        template<typename T>
        using PerAcidData = std::unordered_map<AminoAcid, T>;

        /// Specificity parameters for every amino acid type.
        PerAcidData<SpecificityParams> specificity;

        /**
         * Radii of amino acids (for use in non-full-atomic derivation of a
         * contact map from the native structure).
         */
        PerAcidData<double> radius;

        template<typename T>
        using PerPairData = std::unordered_map<
            std::pair<AminoAcid, AminoAcid>, T>;

        /**
         * A _maximum_ in a distribution of sidechain-sidechain contacts
         * between two given amino acids.
         */
        PerPairData<double> pairwiseMinDist;

        /**
         * An optional MJ matrix of Lennard-Jones potential depths depending on
         * the amino acid types in the pair.
         */
        std::optional<PerPairData<double>> mjMatrix;
    };
}