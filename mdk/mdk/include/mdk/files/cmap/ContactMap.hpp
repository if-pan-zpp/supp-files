#pragma once
#include <vector>

namespace mdk::cmap {
    /**
     * A POD representing a contact map (a.k.a. structured part), as represented
     * in the legacy contact map file): specifically, the location and extent
     * of the structured part, associated native bond and dihedral angles and
     * Go-model contacts.
     *
     * This data structure doesn't provide its own means of parsing
     * (\p cmap::LegacyParser must be used, for example), since a change of
     * format shouldn't influence this POD.
     */
    class ContactMap {
    public:
        /// Length of the structured part.
        int len;

        /// Offset of the structured part (w.r.t. the chain beginning).
        int offset;

        struct Contact {
            /// Residues forming the contact.
            int res[2];

            /// Native length of the contact.
            double dist0;
        };

        /// List of Go-model contacts.
        std::vector<Contact> contacts;

        /**
         * angle[i] denotes the native bond angle value for triple
         * (i-1, i, i+1), relative to the position as indicated by /p off,
         * i.e. if (i+off-1, i+off, i+off+1) relative to the start of the chain.
         * Note: the native angles near the ends of the structured part are
         * undefined.
         */
        std::vector<double> angle;

        /**
         * dihedral[i] denotes the native dihedral angle value for quadruple
         * (i-2, i-1, i, i+1), relative to the position as indicated by /p off,
         * i.e. if (i+off-2, i+off-1, i+off, i+off+1) relative to the start
         * of the chain. Note: the native angles near the ends of the
         * structured part are undefined.
         */
        std::vector<double> dihedral;
    };
}