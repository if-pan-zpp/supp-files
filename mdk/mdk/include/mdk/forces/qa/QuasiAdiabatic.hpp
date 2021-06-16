#pragma once
#include "../NonlocalForce.hpp"
#include "../../model/Model.hpp"
#include "../../files/param/Parameters.hpp"
#include "../../kernels/LennardJones.hpp"
#include "../../kernels/SidechainLJ.hpp"
#include "../../data/Types.hpp"
#include "../../data/Chains.hpp"
#include "../../stats/Stats.hpp"
#include "../../data/Primitives.hpp"
#include <mutex>

namespace mdk {
    /**
     * A struct detailing the present QA contact.
     */
    struct QAContact {
        /// Index of the first residue.
        int i1;

        /// Index of the second residue.
        int i2;

        /// Enum for the status of a QA contact.
        enum class Status: int8_t {
            /// The contact is forming.
            FORMING,
            /// The contact is breaking.
            BREAKING,
            /**
             * The contact has been removed (but since we don't want to
             * actually remove items from the vector, we just mark it this way.)
             */
            REMOVED };

        /// Status of the contact.
        Status status;

        /// Type of contact (BB, BS, SB, SS).
        Stats::Type type;

        /// Time of either formation or when contact has started breaking.
        double t0;
    };

    /**
     * A struct detailing a free pair that can form a QA contact.
     */
    struct QAFreePair {
        /// Index of the first residue.
        int i1;

        /// Index of the second residue.
        int i2;

        /// Enum for the status of the free pair.
        enum class Status: int8_t {
            /// The pair is indeed free.
            FREE,
            /**
             * The pair has formed a QA contact, but since we don't want to
             * outright remove an element from a list, we only mark it as
             * removed.
             */
            TAKEN
        };

        /// Status of the pair.
        Status status;
    };

    /**
     * A "diff" with respect to the QA contact. Since the formation of the
     * contact requires there to be enough "stat slots", and it's possible that
     * in an asynchronous execution more contacts would be formed while they
     * shouldn't, i.e. oversaturate the pool of stat slots, we separate the
     * formation into a "virtual formation" which is done asynchronously and
     * which adds a \p QADiff to a list of diffs, and a synchronous part which
     * adds these diffs one after another, thus preventing an oversaturation.
     */
    struct QADiff {
        /// QA contact to add.
        QAContact cont;

        /// Index of a free pair that now got taken.
        int oldIdx;

        /// Diffs to \p Stat that result from the formation of a QA contact.
        Stat statDiffs[2];

        /**
         * Linear order on \p QADiff structures, done in order to ensure the
         * determinism of the process by sorting the diffs' array.
         */
        bool operator<(QADiff const& other) const {
            return oldIdx < other.oldIdx;
        }
    };

    /**
     * The quasi-adiabatic potential.
     */
    class QuasiAdiabatic: public NonlocalForce {
    public:
        /**
         * Bind the force to the simulation. In particular, we fetch the
         * \p types and \p chains, load sidechain distances etc.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation &simulation) override;

        /**
         * Asynchronous part of the computation. In particular, we compute
         * forces between currently-existing pairs and check which of the
         * free pairs could form contacts.
         * @param dynamics Dynamics object to add potential energy and
         * forces to.
         */
        void asyncPart(Dynamics &dynamics) override;

        /**
         * Synchronous part of the computation. We sort the list of
         * potentially-added contacts and try to add them sequentially.
         * @param dynamics Dynamics object to add potential energy and
         * forces to; here it should be unused unless we want to account for
         * the forces and potential energy of the added contacts.
         */
        void syncPart(Dynamics &dynamics) override;

        /**
         * Action to be performed when the Verlet list is updated. Here we
         * update the lists of pairs in contact and free pairs, preserving the
         * pairs that were in contact in the old list and is present in the
         * new Verlet list.
         */
        void vlUpdateHook() override;

    private:
        /**
         * Generate a spec for the Verlet list.
         * @return Generated spec for the Verlet list.
         */
        vl::Spec spec() const override;

        /**
         * Types of residues, used when selecting the sidechain-sidechain
         * potential.
         */
        Types const *types = nullptr;

        /**
         * Chains data structure, here used mostly to filter the residues
         * that are at the ends of the residues (there $n_i$ and $h_i$ cannot
         * be computed).
         */
        Chains const *chains = nullptr;

        /// Lennard-Jones potential for the backbone-backbone contacts.
        LennardJones bb_lj = LennardJones(5.0 * angstrom, 1.0 * eps);

        /// Lennard-Jones potential for the backbone-sidechain contacts.
        LennardJones bs_lj = LennardJones(6.8 * angstrom, 1.0 * eps);

        /**
         * An array of sidechain L-J potentials, one for every pair of types of
         * amino acids.
         */
        SidechainLJ ss_ljs[AminoAcid::N][AminoAcid::N];

        /// An array of vectors $n_i$, as defined in CPC14.pdf.
        Vectors n;

        /// An array of vectors $h_i$, as defined in CPC14.pdf.
        Vectors h;

        /**
         * A maximum distance (across all types of contacts) between residues
         * for the formation to occur. Has no "physical" sense, is used chiefly
         * for optimization.
         */
        double formationMaxDistSq;

        /**
         * Minimum value of <h_i, r_j> for the formation of
         * a backbone contact (with backbone part for i).
         */
        double hr_abs_min = 0.92;

        /**
         * Minimum value of <h_i, h_j> for the formation of
         * a backbone-backbone contact to occur.
         */
        double hh_abs_min = 0.75;

        /**
         * Maximum value of <n_i, r_j> for the formation of
         * a sidechain contact (with sidechain part for i).
         */
        double nr_max = 0.5;

        /**
         * A scalar factor by which the "usual" minimum formation distance
         * is multiplied.
         */
        double formationTolerance = 1.0;

        /**
         * A timespan during which a formed QA contact "thermalizes", i.e.
         * reaches the full depth of the L-J potential. The increase of the
         * depth is linear.
         */
        double formationTime = 10.0 * tau;

        /**
         * A scalar factor by which the "usual" maximum distance before breaking
         * is multiplied.
         */
        double breakingTolerance = 1.0;

        /**
         * Timespan during which a broken QA contact "dissipates", i.e.
         * loses its force. The decrease in the depth of the L-J potential
         * is linear.
         */
        double breakingTime = 10.0 * tau;

        /**
         * The list of old pairs, with which the new pairs are swapped. This
         * is done in order to not have to allocate new memory each time a
         * Verlet list is regenerated.
         */
        std::vector<QAContact> oldPairs;

        /**
         * A list of current active (or destroyed) QA contacts.
         */
        std::vector<QAContact> pairs;

        /**
         * A list of free pairs, i.e. ones which are not in a QA contact but
         * which potentially could form one.
         */
        std::vector<QAFreePair> freePairs;

        /**
         * A reference to \p Stats variable.
         * Note: this is a non-const reference, and in particular it gets
         * potentially modified in the \p syncPart.
         */
        Stats *stats;

        /**
         * A mutex for the access to the \p qaDiffs list. We have decided not
         * to employ a separate list for every thread since (we estimate) there
         * shouldn't be many conflicts to the access; nevertheless we wish to
         * prevent any concurrent access were one to occur.
         */
        std::mutex qaDiffsMutex;

        /**
         * A list of QA diffs to pass through in the \p syncPart, gets new
         * diffs added to in the formation pass.
         */
        std::vector<QADiff> qaDiffs;

        /**
         * Compute a list of vectors $n_i$ and $h_i$.
         */
        void computeNH();

        /**
         * Perform a geometry check between two residues during the
         * formation pass.
         * @param p Data regarding the pair currently processed.
         * @param diff A \p QADiff structure to modify if the check is
         * successful.
         * @return \p true if the geometry check has succeeded,
         * \p false otherwise.
         */
        bool geometryPhase(vl::PairInfo const& p, QADiff& diff) const;
    };
}
