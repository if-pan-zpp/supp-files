#pragma once
#include "../model/Model.hpp"
#include "../data/Primitives.hpp"

namespace mdk {
    /**
     * This data structure contains all the data pertaining to chains, especially
     * the pairs, triples etc. in the model. It's here due to similarity of
     * computation and some degree of reuse of the same data.
     *
     * The definitions of pairs, triples etc. follow convention from CPC14.pdf.
     */
    class Chains {
    public:
        /**
         * Const ptr to an underlying model (necessarily it must have at least
         * as long a lifetime as this class).
         */
        Model const* model = nullptr;

        /**
         * chainIdx[i] is the index of the chain for residue i, or -1 if a
         * residue is in no chain.
         */
        Integers chainIdx;

        /**
         * isTerminal[i] = 1 if the i'th residue is first or last in the chain
         * (and thus lacks a neighbor), or 0 otherwise.
         */
        Bytes isTerminal;

        /**
         * pairs[i] = 1 if the pair (i, i+1) is entirely in one chain, or
         * 0 otherwise. Result of tuples(2). Used in harmonic tether forcefield.
         */
        Bytes pairs;

        /**
         * triples[i] = 1 if the triple (i-1, i, i+1) is entirely in one chain,
         * or 0 otherwise. Result of tuples(3). Used for QA (in computing n, h)
         * and for bond angle potentials.
         */
        Bytes triples;

        /**
         * quads[i] = 1 if the quadruple (i-2, i-1, i, i+1) is entirely in
         * one chain, or 0 otherwise. Result of tuples(4). Used in dihedral
         * potential.
         */
        Bytes quads;

        /**
         * nativeTriples[i] = 1 if the triple (i-1, i, i+1) is entirely in
         * a structured part, or 0 otherwise. Result of nativeTuples(3). Used
         * in native bond angle potential.
         */
        Bytes nativeTriples;

        /**
         * nativeQuads[i] = 1 if the quadruple (i-2, i-1, i, i+1) is entirely
         * in a structured part, or 0 otherwise. Result of nativeTuples(4).
         * Used in native dihedral potentials.
         */
        Bytes nativeQuads;

        Chains() = default;

        /** Generate a Chains structure for a given model.
         * @param model The model, chain data whereof we wish to obtain.
         */
        explicit Chains(Model const& model);

        /**
         * Check bond distance between residues.
         * @param i1 First residue
         * @param i2 Second residue
         * @param n Bond offset
         * @return true if residues i1 and i2 are either in separate chains or
         * in a single chain but separated by at least n bonds.
         */
        inline bool sepByAtLeastN(int i1, int i2, int n) const {
            return abs(i1 - i2) >= n || (chainIdx[i1] != chainIdx[i2]);
        }

        /**
         * Create an array of tuple indicators.
         * @param k Tuple size
         * @return An array of size model->n, where tuples(k)[i] is 1 if a tuple
         * (i-k+2, ..., i+1) is entirely contained in a chain, and 0 otherwise.
         */
        Bytes tuples(int k) const;

        /**
         * Create an array of native tuple indicators.
         * @param k Tuple size
         * @return An array of size model->n where tuples(k)[i] is 1 if a tuple
         * (i-k+2, ..., i+1) is entirely contained in a structured part, and 0
         * otherwise.
         */
        Bytes nativeTuples(int k) const;
    };
}