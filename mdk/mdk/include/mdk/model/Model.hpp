#pragma once
#include <vector>
#include <string>
#include <optional>
#include <Eigen/Core>
#include "../utils/Units.hpp"
#include "../utils/Random.hpp"
#include "../utils/Topology.hpp"
#include "../utils/ResType.hpp"
#include "../utils/ContactType.hpp"
#include "../files/cmap/ContactMap.hpp"

namespace mdk {
    /**
     * A model object, which serves as one of two objects (other being
     * \p param::Parameters) which instantiate the simulation. Contains
     * residues (with types, messes, positions etc.), chains and structured
     * parts, topological information (i.e. the box shape), along with an array
     * of utilities for the modification of the model, in particular initializing
     * the velocities or preparing the model into a conformation (a line or a
     * self-avoiding walk).
     */
    class Model {
    public:
        /**
         * The number of residues in the model.
         */
        int n = 0;

        struct Residue;
        struct Chain;

        /**
         * A structure representing a single "residue", or a pseudoatom, of the
         * model.
         */
        struct Residue {
            /**
             * Index of the residue in \p residues.
             */
            int idx;

            /**
             * Index of the parent chain, or -1 if it doesn't belong to any
             * chain.
             */
            int chainIdx;

            Vector r, v;
            double mass;
            ResType type;
            /**
             * Native positions of the residues.
             */
            std::optional<Vector> nat_r;
        };
        std::vector<Residue> residues;
        Residue& addResidue(Chain *chain = nullptr);

        struct Chain {
            int idx;
            int start, end;
            /**
             * A list of indices of structured parts that are overlayed over
             * the chain.
             */
            std::vector<int> structuredParts;
        };
        std::vector<Chain> chains;
        Chain& addChain();

        struct Contact {
            int idx;
            int res[2];
            ContactType type;
            double dist0;
        };
        std::vector<Contact> contacts;
        Contact& addContact();

        struct StructuredPart {
            int idx;
            int off, len;
            std::vector<double> angle, dihedral;
        };
        std::vector<StructuredPart> structuredParts;
        StructuredPart& addSP();

        /**
         * Topology, or the box size.
         */
        Topology top;

    public:
        void morphIntoLine();

        /**
         * Morphing the model into a self-avoiding walk. A custom version
         * with clarity in mind, doesn't replicate the positions from the
         * Fortran code which makes testing difficult.
         * @param rand Random number source to use.
         * @param useTop Whether to create and use the topology during the
         * creation process.
         * @param density Target density of the box. 0 means an infinite box?
         * @param minDist Minimum distance for the conformation to be considered
         * intersecting, and therefore invalid.
         *
         * Note: it doesn't use native bond angles correctly at the moment.
         */
        void morphIntoSAW(Random& rand,
            bool useTop = false,
            double density = 1e-4 * atom / pow(angstrom, 3.0),
            double minDist = 4.56 * angstrom);

        /**
         * Morphing the model into a self-avoiding walk. A version which
         * replicates the behavior of the Fortran code but in a more clear
         * way.
         * @param rand Random number source to use.
         * @param useTop Whether to create and use the topology during the
         * creation process.
         * @param density Target density of the box. 0 means an infinite box?
         * @param minDist Minimum distance for the conformation to be considered
         * intersecting, and therefore invalid.
         * @param nativeBondLen Whether to derive the bond lengths from the
         * native structure. By default the bond lengths are 3.8 angstrem; if
         * \p nativeBondLen is true, the bond length is taken to be a mean bond
         * length in the native structure.
         */
        void legacyMorphIntoSAW(Random& rand,
            bool useTop = false,
            double density = 1e-4 * atom / pow(angstrom, 3.0),
            double minDist = 4.56 * angstrom,
            bool nativeBondLen = false);

        /**
         * Morphing the model into a self-avoiding walk. A version which
         * replicates (one-by-one) the code of the Fortran version. It is also
         * less intuitive because of it.
         * @param rand Random number source to use.
         * @param useTop Whether to create and use the topology during the
         * creation process.
         * @param density Target density of the box. 0 means an infinite box?
         * @param minDist Minimum distance for the conformation to be considered
         * intersecting, and therefore invalid. Note, from what I see the
         * intersection is checked only for each chain internally, and not
         * across the chains.
         * @param bond Default bond length to use.
         */
        void exactLegacySAW(Random& rand,
            bool useTop = false,
            double density = 0.0 * atom / pow(angstrom, 3.0),
            double minDist = 4.56 * angstrom,
            double bond = 3.8 * angstrom);

        /**
         * Initialize velocities of the residues, assuming they are placed in
         * a medium with a given temperature.
         * @param rand Random number source to use.
         * @param temperature Temperature of the medium.
         * @param useMass Whether to account for different masses of the
         * residues.
         */
        void initVelocity(Random& rand,
                          double temperature,
                          bool useMass = true);

        StructuredPart& addContactMap(cmap::ContactMap const& contactMap);
        void addCMapContacts(cmap::ContactMap const& contactMap, Chain& chain);

        /**
         * Set the masses of all residues to be equal to the mean of the
         * current masses.
         */
        void useAverageMasses();

    private:
        /**
         * Generates the pairs of residues separated by at least two bonds,
         * used in the SAW procedures.
         * @return Pairs of residues.
         */
        std::vector<std::pair<Residue*, Residue*>> nonlocalPairs();
    };
}
