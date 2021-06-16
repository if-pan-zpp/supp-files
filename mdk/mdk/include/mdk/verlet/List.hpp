#pragma once
#include "../system/State.hpp"
#include "../utils/Units.hpp"
#include "../simul/SimulVar.hpp"
#include "../data/Chains.hpp"
#include "Spec.hpp"

namespace mdk {
    class NonlocalForce;
}

namespace mdk::vl {
    /**
      * A struct containing data about pairs; this is stored in this fashion
       * chiefly to avoid having to pass all of these as arguments.
       */
    struct PairInfo {
        /**
         * Index of the first residue.
         */
        int i1;

        /**
         * Index of the second residue.
         */
        int i2;

        /**
         * Normalized vector connecting the two residues. May be non-trivial,
         * i.e. not just r_{i_1} - r_{i_2} but the PBC-corrected version thereof.
         */
        Vector unit;

        /**
         * Norm of the vector connecting the two residues.
         */
        double norm;
    };

    /**
     * A Verlet list.
     */
    class List: public SimulVar {
    private:
        State const *state = nullptr;
        Chains const* chains = nullptr;

        /**
         * Time from when the list was last reconstructed.
         */
        double t0 = 0.0;

        /**
         * Positions of the residues from when the list was last
         * reconstructed.
         */
        Vectors r0;

        /**
         * Box shape from when the list was last reconstructed.
         */
        Topology top0;

        /**
         * A list of nonlocal forces, saved in order to invoke their
         * update hooks when the list is reconstructed.
         */
        std::vector<NonlocalForce*> forces;

        bool initial = false;

        /**
         * A maximal cutoff among registered nonlocal forces.
         */
        double cutoff = 0.0 * angstrom;

        /**
         * An extra "buffer". The list by default contains elements that may be
         * outside the maximal cutoff, which allows us to not have to
         * reconstruct the verlet list at every time step, at the cost of
         * spurious distance comparisons.
         */
        double pad = 10.0 * angstrom;

        /**
         * Minimal bond-distance between the residues among registered nonlocal
         * forces.
         */
        int minBondSep = 0;

        bool needToReset() const;

        /**
         * Updates the Verlet list in the "legacy fashion", i.e. going through
         * a list of pairs and checking the distances.
         */
        void update();

        /**
         * Dimensions of the grid of cells used in the cell version of the
         * computation.
         */
        Eigen::Vector3i grid;

        /**
         * A list of size equal to the number of cells in the grid, where a
         * nonnegative number indicates an index of the residue which is first
         * in the linked list of residues in a given cell. If the cell is empty,
         * the value is -1.
         */
        std::vector<int> first;

        /**
         * A list of size equal to the number of cells in the grid, where a
         * nonnegative number indicates an index of the residue which is last
         * in the linked list of residues in a given cell. If the cell is empty,
         * the value is -1.
         */
        std::vector<int> last;

        /**
         * A list of size equal to the number of residues, where a nonnegative
         * number indicates an index of a residue which is next in the linked
         * list associated with the cell in which the residue is placed.
         */
        std::vector<int> next;

        double effCutoff, effCutoffSq;

        /**
         * Computes the flattened index of a cell in the grid.
         * @param loc Unflattened index.
         * @return Flattened index.
         */
        int indexOf(Eigen::Vector3i const& loc);

        void perPair(int c1, int c2);
        void perCell(int c1);

        /**
         * Updates the Verlet list, but instead of going through each pair we
         * first place the residues into cells of (rough) size \p effCutoff
         * and then only take the pairs from the neighboring cells, which
         * reduces the computational cost from O(N^2) to O(N) (although the
         * constant behind O(N) may be significant). In practice this version
         * is much faster, especially for large numbers of residues.
         */
        void updateGrid();

    public:
        /**
         * Register a nonlocal force, in particular adjust the current specs
         * so as to accomodate the newly added force (for example increase the
         * cutoff distance), and add it to the internal list of the forces whose
         * hooks to call when a list is reconstructed.
         * @param force Nonlocal force to register.
         * @param spec Spec with which to register the force.
         */
        void registerNF(NonlocalForce& force, Spec const& spec);

        /**
         * The list of pairs. The pairs are ordered (if a pair [i, j] is in
         * this list, then i < j).
         */
        Pairs pairs;

        void bind(Simulation& simulation) override;

        /**
         * Check whether the list needs updating, and if it does update it
         * and invoke the relevant update hooks for the non-local forces.
         */
        void check();
    };
}
