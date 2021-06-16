#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <optional>
#include "../../model/Model.hpp"
#include "../param/Parameters.hpp"

namespace mdk::pdb {
    /**
     * PDB model object. As opposed to \p Model, here we store atomic data,
     * chiefly for the use in the derivation of full-atomic contact maps. Also,
     * non-full-atomic contact maps may be derived; \p Model doesn't contain such
     * a method because it's ill-defined for models obtained from \p Sequence
     * file, which doesn't have native residue positions.
     */
    class Model {
    public:
        struct Atom;
        struct Residue;
        struct Chain;

        /// A struct containing data for a single atom.
        struct Atom {
            /// PDB serial number of an atom. Zero-indexed.
            int serial;

            /**
             * Index of an atom in the \p Residue::atoms list of the
             * parent residue, or -1 if it's not a part of any residue.
             */
            int idxInRes;

            /// Parent residue.
            Residue *res;

            /// (Native) position of the atom.
            Eigen::Vector3d r;

            /// Type of an atom (in the \p std::string form).
            std::string type;
        };

        /**
         * A map from the serial number of an atom to the \p Atom structure.
         * Note that the order of serial numbers is (1) not necessarily
         * consecutive, (2) not starting from 0, in some degenerate cases.
         * Also, adding new atoms doesn't invalidate pointers or references.
         */
        std::unordered_map<int, Atom> atoms;

        /**
         * Add an atom to a residue.
         * @param serial Serial number of an atom
         * @param res Parent residue, or \p nullptr
         * @return Reference to the added atom.
         */
        Atom& addAtom(int serial, Residue *res = nullptr);

        /// A struct containing data for a PDB residue.
        struct Residue {
            /// Serial number of a residue. 0-indexed.
            int serial;

            /**
             * Index in the \p Chain::residues of the parent chain, or -1 if
             * it's not in any chain.
             */
            int idxInChain;

            /// Parent chain.
            Chain *chain;

            /// Type of residue (in \p std::string form)
            std::string type;

            /// List of constituent atoms for the residue.
            std::vector<Atom*> atoms;

            /**
             * Find an atom belonging to the residue by the type. (Used
             * primarily for finding CG atoms in SSBONDs).
             * @param type Atom type to find.
             * @return Pointer to found atom, or \p nullptr if such an atom
             * has not been found.
             */
            Atom* find(std::string const& type);
        };

        /**
         * A map from residue serial number to Residue struct. Note that the
         * order of serial numbers is (1) not necessarily consecutive,
         * (2) not starting from 0, in some degenerate cases. Also, adding
         * new residues doesn't invalidate pointers or references.
         */
        std::unordered_map<int, Residue> residues;

        /**
         * Add a residue to a chain.
         * @param serial Serial number of a residue to add.
         * @param chain Parent chain or \p nullptr
         * @return Reference to the added residue.
         */
        Residue& addResidue(int serial, Chain *chain = nullptr);

        /// A struct containing PDB chain data.
        struct Chain {
            /// Serial (alphanumeric) character of a chain.
            char serial;

            /// Constituent resides of the chain. Ordered.
            std::vector<Residue*> residues;
        };

        /// A map from the chain serial character to the chain.
        std::unordered_map<char, Chain> chains;

        /**
         * Add a chain with a specified serial character.
         * @param serial Serial character of the chain.
         * @return Reference to the added chain.
         */
        Chain& addChain(char serial);

        /**
         * A struct containing contacts. Technically this isn't represented
         * 1-to-1 in the PDB data, but an unification of SSBONDs and LINKs is
         * useful, especially for the reduction of the \p pdb::Model into \p
         * Model.
         */
        struct Contact {
            /// Index of the contact in \p contacts vector.
            int idx;

            /// Pointers to atoms that comprise the contact.
            Atom* atom[2];

            /** Type of contact. Takes (at this moment) values "SSBOND",
             * "NAT", "NAT_BB", "NAT_BS", "NAT_SB", "NAT_SS".
             */
            std::string type;

            /// Native contact length.
            double dist0;
        };

        /**
         * A vector of contacts. Note that references _are_ (potentially)
         * invalidated when adding new contacts.
         */
        std::vector<Contact> contacts;

        /**
         * Add a new contact.
         * @return Reference to the added contact.
         */
        Contact& addContact();

        /// "Topology", i.e. box shape for PBC.
        Topology top;

    public:
        Model() = default;

        /**
         * Construct a \p pdb::Model from a \p Model. Doesn't extrapolate extra
         * atoms, i.e. only CA atoms are in the final model.
         * @param coarse \p Model from which to construct the PDB model.
         */
        explicit Model(mdk::Model const& coarse);

        /**
         * Adds contacts based on the overlap of atoms. This method also
         * derives types of contacts based on whether the atoms are in the
         * backbone or the sidechain of an amino acid.
         */
        void addContactsFromAtomOverlap();

        /**
         * Adds contacts based on approximate residue overlap, given the
         * \p aminoAcidRadii data from the parameter file.
         * @param params Parameter file with amino acid radii data.
         */
        void addContactsFromResOverlap(param::Parameters const& params);

        /**
         * Construct a \p mdk::Model from the PDB model. Collapses residues to
         * a single \p Model residue and copies the contacts (if two atoms are
         * connected, the residues containing them are connected in the final
         * model).
         * @return A reduced ("coarsened") model.
         */
        mdk::Model coarsen();
    };
}
